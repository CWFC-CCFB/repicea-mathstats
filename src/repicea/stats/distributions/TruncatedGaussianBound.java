package repicea.stats.distributions;

import java.io.Serializable;

import repicea.math.Matrix;
import repicea.stats.distributions.utility.GaussianUtility;

@SuppressWarnings("serial")
public class TruncatedGaussianBound extends BasicBound implements Serializable {

	/**
	 * Ensure the instance provide the information for properly set the Truncated GaussianBound
	 * @author Mathieu Fortin - July 2022
	 */
	public interface TruncatedGaussianBoundCompatible {
		/**
		 * Return the mean of the original distribution (before truncating)
		 * @return a double
		 */
		public double getMuValue();
		
		/**
		 * Return the variance of the original distribution (before truncating)
		 * @return
		 */
		public double getSigma2Value();
	}
	
	
	private boolean isCompletelySet;
	private double pdfValueOnStandardNormal;
	private double cdfValue;
	private double standardizedValue;
	private Object lock = new Object();

	private final TruncatedGaussianBoundCompatible owner;
	
	public TruncatedGaussianBound(TruncatedGaussianBoundCompatible owner, boolean isUpperBound) {
		super(isUpperBound);
		this.owner = owner;
	}

	@Override
	public void setBoundValue(Matrix value) {
		super.setBoundValue(value);
		isCompletelySet = false;
	}
	
	public double getPdfValueOnStandardNormal() {
		synchronized(lock) {
			if (!isCompletelySet) {
				update();
			}
		}
		return pdfValueOnStandardNormal;
	}
	
	public double getCdfValue() {
		synchronized(lock) {
			if (!isCompletelySet) {
				update();
			} 
		}
		return cdfValue;
	}
	
	private double getBoundValueDouble() {
		return getBoundValue().getValueAt(0, 0);
	}

	private void update() {
		if (getBoundValue() == null) {
			pdfValueOnStandardNormal = 0d;
			standardizedValue = Double.POSITIVE_INFINITY;
			if (isUpperBound()) {
				cdfValue = 1d;
			} else {
				cdfValue = 0d;
			}
		} else {
			double std = Math.sqrt(owner.getSigma2Value());
			standardizedValue = (getBoundValueDouble() - owner.getMuValue()) / std;
			pdfValueOnStandardNormal = GaussianUtility.getProbabilityDensity(standardizedValue);
			cdfValue = GaussianUtility.getCumulativeProbability(standardizedValue);
		}
		isCompletelySet = true;
	}
	
	public double getStandardizedValue() {
		synchronized(lock) {
			if (!isCompletelySet) {
				update();
			} 
		}
		return standardizedValue;
	}
	
	public void reset() {
		isCompletelySet = false;
	}
}