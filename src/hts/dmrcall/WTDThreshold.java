package hts.dmrcall;

import java.util.HashMap;
import java.util.Map;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

import cern.jet.stat.Gamma;

import utils.AbstractLineProcessor;

@App(overview="Apply p-value thresholds to WindowedTagDiff data files", generateStub=true)
public class WTDThreshold extends AbstractLineProcessor {
	private double rat = 0.5;
	private double threshold = 0.001;
	private int validated = 0;
	private int reported = 0;
	private Map<Integer,double[]> pCache = new HashMap<Integer, double[]>();
	private boolean printValidated = false;
	private int stupidDepth = 100000000;
	
	@Option(help="...", optional=true)
	public void setStupidDepth(int i) {
		this.stupidDepth = i;
	}
	
	@Option(help="...", optional=true)
	public void setPrintValidated(boolean b) {
		this.printValidated = b;
	}
	
	@Option(help="...", optional=true)
	public void setThreshold(double d) {
		this.threshold = d;
	}
	
	@Option(help="...", optional=true)
	public void setRat(double d) {
		this.rat = d;
	}
	
	private double[] getPvals(int N)
	{
		double[] pv = pCache.get(N);
		if (pv == null) {
			pv = new double[N + 1];
			for (int n = 0; n <= N; ++n) {
				double p = 0;
				double mid = rat * N;
				if (n > mid) {
					for (int x = n; x <= N; ++x) {
						p += binom(rat, x, N);
					}
				} else {
					for (int x = n; x >= 0; --x) {
						p += binom(rat, x, N);
					}
				}
				pv[n] = p;
			}
			pCache.put(N, pv);
		}
		return pv;
	}
	
	private boolean isValid(double[] pv) {
		for (double d : pv) {
			if (d < threshold) {
				return true;
			}
		}
		return false;
	}
	
	public void processTokens(String[] toks) {
		int m = Integer.parseInt(toks[0]);
		int n = Integer.parseInt(toks[1]);
		int N = n + m;
		
		if (N > stupidDepth) {
			return;
		}
		
		double[] pv = getPvals(N);
		if (isValid(pv)) {
			++validated;
		}
		
		if (pv[n] < threshold || (printValidated && isValid(pv))) {
			++reported;
			for (String tok : toks) {
				System.out.print(tok);  System.out.print('\t');
			}
			System.out.printf("%g\t%g\t%s%n", pv[n], -Math.log10(pv[n]), n < (rat * N) ? "hypo" : "hyper");
		}
	}
	
	public void post() {
		System.err.printf("Estimated FDR: %g%n", (threshold * validated) / reported);
	}
	
    private static double binom(double prob, int n, int N) {
        try {
            return Math.exp(
                Gamma.logGamma(N + 1) -
                Gamma.logGamma(n + 1) - 
                Gamma.logGamma(N - n + 1) +
                n * Math.log(prob) +
                (N - n) * Math.log(1 - prob)
            );
        } catch (ArithmeticException ex) {
            System.err.println("Error calculating binom(" + prob + ", " + n + ", " + N + ")");
            throw ex;
        }
    }
}
