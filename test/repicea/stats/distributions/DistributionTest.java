package repicea.stats.distributions;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimates.MonteCarloEstimate;

public class DistributionTest {

	/*
	 * This test uses a Multivariate Gaussian distribution as a reference and compute a Monte Carlo estimator. The mean and the variances are then compared.
	 */
	@Test
	public void stochasticSimulationFromGaussianDistribution() {
		int nbReal = 1000000;
		EmpiricalDistribution npDist = new EmpiricalDistribution(); 
		Matrix mean = new Matrix(2,1);
		mean.setValueAt(0, 0, 2d);
		mean.setValueAt(1, 0, 3d);
		SymmetricMatrix variance = new SymmetricMatrix(2);
		variance.setValueAt(0, 0, 0.5);
		variance.setValueAt(1, 0, 0.4);
		variance.setValueAt(1, 1, 1d);
		
		GaussianEstimate estimate = new GaussianEstimate(mean, variance);
		
		for (int i = 0; i < nbReal; i++) {
			npDist.addRealization(estimate.getRandomDeviate());
		}
		
		Matrix simulatedMean = npDist.getMean();
		Matrix res = simulatedMean.subtract(mean);
		double sse = res.transpose().multiply(res).getValueAt(0, 0);
		System.out.println("Squared difference of the means = " + sse);
		
		Assert.assertEquals("Testing the means", 0d, sse, 1E-4);

		Matrix simulatedVariances = npDist.getVariance();
		Matrix diff = simulatedVariances.subtract(variance);
		sse = diff.elementWiseMultiply(diff).getSumOfElements();
		System.out.println("Squared difference of the variances = " + sse);

		Assert.assertEquals("Testing the variances", 0d, sse, 1E-4);

	}
	
	/*
	 * This test uses a univariate truncated Gaussian distribution as a reference and compute a Monte Carlo estimator. The mean and the variances are then compared.
	 */
	@Test
	public void stochasticSimulationFromTruncatedStandardGaussianDistribution() {
		int nbReal = 1000000;
		
		TruncatedGaussianDistribution distribution = new TruncatedGaussianDistribution();
		distribution.setUpperBoundValue(new Matrix(1,1));
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		
		for (int i = 0; i < nbReal; i++) {
			estimate.addRealization(distribution.getRandomRealization());
		}
		
		Matrix simulatedMean = estimate.getMean();
		Matrix simulatedVariance = estimate.getVariance();
		Matrix expectedMean = distribution.getMean();
		Matrix expectedVariance = distribution.getVariance();
		Assert.assertEquals("Testing the means", expectedMean.getValueAt(0, 0), simulatedMean.getValueAt(0, 0), 2E-3);
		Assert.assertEquals("Testing the variances", expectedVariance.getValueAt(0, 0), simulatedVariance.getValueAt(0, 0), 3E-3);
	}

	
	/*
	 * This test uses a univariate truncated Gaussian distribution as a reference and compute a Monte Carlo estimator. The mean and the variances are then compared.
	 */
	@Test
	public void stochasticSimulationFromTruncatedGaussianDistribution() {
		int nbReal = 1000000;
		
		TruncatedGaussianDistribution distribution = new TruncatedGaussianDistribution(10, 20);
		Matrix upperBound = new Matrix(1,1);
		upperBound.setValueAt(0, 0, 12);
		distribution.setUpperBoundValue(upperBound);
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		
		for (int i = 0; i < nbReal; i++) {
			estimate.addRealization(distribution.getRandomRealization());
		}
		
		double simulatedMean = estimate.getMean().getValueAt(0, 0);
		double simulatedVariance = estimate.getVariance().getValueAt(0, 0);
		double expectedMean = distribution.getMean().getValueAt(0, 0);
		double expectedVariance = distribution.getVariance().getValueAt(0, 0);
		System.out.println("Expected mean = " + expectedMean + ", actual mean = " + simulatedMean);
		Assert.assertEquals("Testing the means", expectedMean, simulatedMean, 1.5E-2);
		System.out.println("Expected variance = " + expectedVariance + ", actual variance = " + simulatedVariance);
		Assert.assertEquals("Testing the variances", expectedVariance, simulatedVariance, 5E-2);
	}

}
