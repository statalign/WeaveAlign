package test;

import java.util.Arrays;
import java.util.Locale;

import mpd.MpdMain;

public class NumOpt {
	
	boolean autoLogRange = true;
	
	boolean[] logRange;
	
	double[] from, to;
	double[][] curRange;
	
	double[] curPars;
	double[] results;
	
	
	public double[] optimise(NumOptTask task, int steps) {
		int ranl = task.getFrom().length == 1 ? 5 : 4;	// range length
		return optimise(task, steps, ranl);
	}
	
	public double[] optimise(NumOptTask task, int steps, int ranl) {
		from = Arrays.copyOf(task.from, task.from.length);
		to = Arrays.copyOf(task.to, task.to.length);
		int pars = from.length;
		
		for(int i = 0; i < steps; i++) {
			System.out.println("Step "+(i+1));
			
			// generate ranges
			for(int j = 0; j < pars; j++) {
				if(autoLogRange && from[j] > 0 && to[j]/from[j] > 4)
					curRange[j] = logRange(from[j], to[j], ranl);
				else
					curRange[j] = range(from[j], to[j], ranl);
			}
			runCurRange(0);
		}
		return null;
	}
	
	private void runCurRange(int par) {
		if(par == from.length) {
			for(int i = 0; i < curRange[par].length; i++)
				;
		}
	}

	private static double[] range(double from, double to, int len) {
		double[] r = new double[len];
		double step = (to-from)/(len-1);
		for(int i = 0; i < len; i++)
			r[i] = from+i*step;
		return r;
	}

	private static double[] logRange(double from, double to, int len) {
		double[] r = new double[len];
		double lfrom = Math.log(from), lto = Math.log(to);
		double step = (lto-lfrom)/(len-1);
		for(int i = 0; i < len; i++)
			r[i] = Math.exp(lfrom+i*step);
		return r;
	}

	public static void main(String[] args) {
		double[] r = logRange(0.62, 0.78, 8);
		double[] res = new double[r.length];
		for(int i = 0; i < r.length; i++) {
			String arg = String.format(Locale.ENGLISH, "dm2run6.log -mod dm2run6.noncons.mod -rho=0.267 -hmm 2.hmm -g=0 -pseq=droSim1 -gr=%f,2.6",r[i]);
			MpdMain.main(arg.split(" "));
			res[i] = MpdMain.lastDataProb;
		}
		for(int i = 0; i < r.length; i++) {
			System.out.println(String.format(Locale.ENGLISH, "%.3f\t%.3f", r[i], res[i]));
		}
	}
	
}