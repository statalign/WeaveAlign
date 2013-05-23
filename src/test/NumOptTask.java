package test;

public abstract class NumOptTask {
	
	double[] from, to;		// parameter optimisation ranges
	
	abstract double eval(double[] pars);
	
	public double[] getFrom() {
		return from;
	}
	
	public double[] getTo() {
		return to;
	}
}
