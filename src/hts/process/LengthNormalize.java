package data.hts.process;

import gff.GFFUtils;

import io.IOTools;

import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFRecord;
import org.bjv2.util.cli.App;

import utils.Collects;

@App(overview="Match the length distributions between two GFF files", generateStub=true)
public class LengthNormalize {
	private static class Record {
		public final String chr;
		public final int min;
		public final int max;
		
		public Record(String chr, int min, int max) {
			this.chr = chr;
			this.min = min; this.max = max;
		}
	}

	private Map<Integer,List<Record>> load(String n)
		throws Exception
	{
		final Map<Integer,List<Record>> m = new TreeMap<Integer, List<Record>>();
		GFFUtils.parser().parse(IOTools.nameBufferedReader(n), new GFFDocumentHandler() {
			public void commentLine(String comment) {
			}

			public void endDocument() {
			}

			public void recordLine(GFFRecord record) {
				Record r = new Record(record.getSeqName(), record.getStart(), record.getEnd());
				int l = r.max - r.min + 1;
				Collects.pushOntoMap(m, l, r);
			}

			public void startDocument(String locator) {
			}
		});
		return m;
	}
	
	public void main(String[] args) 
		throws Exception
	{
	    if (args.length != 4) {
		System.err.println("Usage java hts.process.LengthNormalizeApplication in1.gff in2.gff out1.gff out2.gff");
		return;
	    }

		Map<Integer,List<Record>> m1 = load(args[0]);
		Map<Integer,List<Record>> m2 = load(args[1]);
		
		PrintWriter pw1 = new PrintWriter(new FileWriter(args[2]));
		PrintWriter pw2 = new PrintWriter(new FileWriter(args[3]));
		for (Integer l : m1.keySet()) {
			if (!m2.containsKey(l)) {
				continue;
			}
			List<Record> l1 = m1.get(l);
			List<Record> l2 = m2.get(l);
			int ml = Math.min(l1.size(), l2.size());
			Collections.shuffle(l1); Collections.shuffle(l2);
			write(pw1, l1, ml);
			write(pw2, l2, ml);
		}
		pw1.close(); pw2.close();
	}
	
	private void write(PrintWriter pw, List<Record> l, int ml)
	{
		for (int i = 0; i < ml; ++i) {
			Record r = l.get(i);
			pw.printf("%s\tln\tread\t%d\t%d\t.\t.\t.%n", r.chr, r.min, r.max);
		}
	}
}
