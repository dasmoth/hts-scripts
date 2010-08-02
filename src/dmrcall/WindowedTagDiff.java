package dmrcall;

import io.IOTools;

import java.io.File;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import gfftools.AbstractTransducer;
import gfftools.GFFUtils;

import net.derkholm.nmica.maths.MathsTools;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.SimpleGFFRecord;
import org.bjv2.util.SmallMap;
import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import cern.jet.stat.Gamma;

@App(overview="Count tag-depth differences between two samples", generateStub=true)
public class WindowedTagDiff {
	private int rez = 1000;
	private int step = 1000;
	private int min = 0;
	private int max = Integer.MIN_VALUE;
	private String chr = null;
	
	private int[] gffToTagArray(File f) 
		throws Exception
	{
		final int[] counts = new int[50000000];
		GFFUtils.parser().parse(IOTools.fileBufferedReader(f), new GFFDocumentHandler() {
			public void commentLine(String comment) {
			}

			public void endDocument() {
			}

			public void recordLine(GFFRecord r) {
				if (chr != null && !chr.equals(r.getSeqName())) {
					return;
				}
				
				int rmin = r.getStart();
				int rmax = r.getEnd();
				
				int centerb = ((rmin + rmax) / 2) / step;
				int minb = centerb;
				while (((minb - 1) * step + rez) >= rmin) {
					minb--;
				}
				int maxb = centerb;
				while (((maxb + 1) * step + 1) <= rmax) {
					maxb++;
				}
				
				minb = Math.max(0, minb); maxb = Math.max(0, maxb);
				
				min = Math.min(min, minb);
				max = Math.max(max, maxb);
				for (int b = minb; b <= maxb; ++b) {
					++counts[b];
				}
			}

			public void startDocument(String locator) {
			}
			
		});
		return counts;
	}

	@Option(help="...", optional=true)
	public void setChr(String s) {
		this.chr = s;
	}
	
	
	@Option(help="...", optional=true)
	public void setRez(int i) {
		this.rez = i;
	}
	
	@Option(help="...", optional=true)
	public void setStep(int i) {
		this.step = i;
	}
	
	
	public void main(String[] args)
		throws Exception
	{
		int[] counts0 = gffToTagArray(new File(args[0]));
		int[] counts1 = gffToTagArray(new File(args[1]));
		SimpleGFFRecord r = new SimpleGFFRecord();
		
		for (int b = min; b <= max; ++b) {
			int c0 = counts0[b];
			int c1 = counts1[b];
			
			System.out.printf("%d\t%d\t%d\t%d\t%g\t%g\t%s\t%d\t%d%n", c0, c1, c0 + c1, c1 - c0, Math.random() + (c0 + c1), Math.random() + (c1 - c0) - 0.5, chr, b * step + 1, b * step + rez);
		}
		
		int total0 = 0, total1 = 0;
		for (int b = min; b <= max; ++b) {
			total0 += counts0[b];
			total1 += counts1[b];
		}
		
		System.err.printf("Expected rat: %g%n", (1.0 * total1) / (total1 + total0));
	}

}
