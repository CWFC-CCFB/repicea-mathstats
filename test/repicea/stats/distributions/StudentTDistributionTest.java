package repicea.stats.distributions;


import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.estimates.ConfidenceInterval;
import repicea.stats.estimates.MonteCarloEstimate;

public class StudentTDistributionTest {

	
	@Test
	public void testRandomNumberGenerationWith3DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(3);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().getValueAt(0, 0);
		double variance = estimate.getVariance().getValueAt(0, 0);
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().getValueAt(0, 0);
		double quantile975 = ci.getUpperLimit().getValueAt(0, 0);
		double expectedMean = dist.getMean().getValueAt(0, 0);
		double expectedVariance = dist.getVariance().getValueAt(0, 0);
		System.out.println("Expected mean " + expectedMean + "; Actual mean " + mean);
		Assert.assertEquals("Testing the mean", expectedMean, mean, 0.003);
		System.out.println("Expected variance " + expectedVariance + "; Actual variance " + variance);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 0.8);
		double expectedQuantile025 = -3.182446;
		System.out.println("Expected quantile 025 " + expectedQuantile025 + "; Actual quantile 025 " + quantile025);
		Assert.assertEquals("Testing quantile 0.025", expectedQuantile025, quantile025, 5E-2);
		double expectedQuantile975 = 3.182446;
		System.out.println("Expected quantile 975 " + expectedQuantile975 + "; Actual quantile 975 " + quantile975);
		Assert.assertEquals("Testing quantile 0.975", expectedQuantile975, quantile975, 5E-2);
	}

	@Test
	public void testRandomNumberGenerationWith10DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(10);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().getValueAt(0, 0);
		double variance = estimate.getVariance().getValueAt(0, 0);
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().getValueAt(0, 0);
		double quantile975 = ci.getUpperLimit().getValueAt(0, 0);
		double expectedMean = dist.getMean().getValueAt(0, 0);
		double expectedVariance = dist.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing the mean", expectedMean, mean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, -2.228139, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, 2.228139, 5E-2);
	}

	@Test
	public void testRandomNumberGenerationWith20DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		StudentTDistribution dist = new StudentTDistribution(20);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double mean = estimate.getMean().getValueAt(0, 0);
		double variance = estimate.getVariance().getValueAt(0, 0);
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().getValueAt(0, 0);
		double quantile975 = ci.getUpperLimit().getValueAt(0, 0);
		double expectedMean = dist.getMean().getValueAt(0, 0);
		double expectedVariance = dist.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing the mean", expectedMean, mean, 2E-3);
		Assert.assertEquals("Testing the variance", expectedVariance, variance, 5E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, -2.085963, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, 2.085963, 5E-2);
	}
	

	@Test
	public void testRandomNumberGenerationWithNonCentered10DegreesOfFreedom() {
		
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Matrix mean = new Matrix(1,1);
		mean.setValueAt(0, 0, 10d);
		SymmetricMatrix variance = new SymmetricMatrix(1);
		variance.setValueAt(0, 0, 20d);
		StudentTDistribution dist = new StudentTDistribution(mean, variance, 10);
		for (int i = 0; i < 5000000; i++) {
			estimate.addRealization(dist.getRandomRealization());
		}
		double actualMean = estimate.getMean().getValueAt(0, 0);
		double actualVariance = estimate.getVariance().getValueAt(0, 0);
		ConfidenceInterval ci = estimate.getConfidenceIntervalBounds(.95);
		double quantile025 = ci.getLowerLimit().getValueAt(0, 0);
		double quantile975 = ci.getUpperLimit().getValueAt(0, 0);
		double expectedMean = dist.getMean().getValueAt(0, 0);
		double expectedVariance = dist.getVariance().getValueAt(0, 0);
		double expectedQuantile025 = -2.228139 * Math.sqrt(20) + 10;
		double expectedQuantile975 = 2.228139 * Math.sqrt(20) + 10;
		Assert.assertEquals("Testing the mean", expectedMean, actualMean, 1E-2);
		Assert.assertEquals("Testing the variance", expectedVariance, actualVariance, 6E-2);
		Assert.assertEquals("Testing quantile 0.025", quantile025, expectedQuantile025, 5E-2);
		Assert.assertEquals("Testing quantile 0.975", quantile975, expectedQuantile975, 5E-2);
	}


}
