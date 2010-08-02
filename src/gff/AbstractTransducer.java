package gff;

import io.IOTools;

import java.io.BufferedReader;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;

import org.biojava.bio.program.gff.GFFDocumentHandler;
import org.biojava.bio.program.gff.GFFParser;
import org.biojava.bio.program.gff.GFFRecord;
import org.biojava.bio.program.gff.GFFWriter;

public class AbstractTransducer implements GFFDocumentHandler {
	private GFFWriter gffw;
	
	/**
	 * @param args
	 */
	public void main(String[] args) 
		throws Exception
	{
		BufferedReader br = IOTools.inputBufferedReader(args);
		
		PrintWriter pw = new PrintWriter(new OutputStreamWriter(System.out));
		gffw = new GFFWriter(pw);
		GFFUtils.parser().parse(br, this);
		pw.flush();
	}

	public void commentLine(String comment) {
		gffw.commentLine(comment);
	}

	public void endDocument() {
		gffw.endDocument();
	}

	public void recordLine(GFFRecord record) {
		gffw.recordLine(record);
	}

	public void startDocument(String locator) {
		gffw.startDocument(locator);
	}

}
