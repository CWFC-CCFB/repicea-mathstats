package repicea.stats.estimates;

import repicea.stats.Distribution;


public interface DistributionProvider<D extends Distribution<?,?>>{
	
	/**
	 * Provide the assumed distribution for the random variable.
	 * @return a Distribution-derived instance
	 */
	public D getDistribution();
	
}
