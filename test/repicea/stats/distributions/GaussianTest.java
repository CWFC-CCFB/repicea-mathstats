package repicea.stats.distributions;

import static org.junit.Assert.assertEquals;

import java.security.InvalidParameterException;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.utility.GaussianUtility;

public class GaussianTest {

	@Test
	public void bivariateCumulativeProbabilitiesTest() {
		double x1 = -0.2876339;
		double x2 = 1.125041;
		double rho = -0.6018568;
		double biv11 = GaussianUtility.getBivariateCumulativeProbability(x1, x2, false, false, rho);
		assertEquals(0.2832904938411274, biv11, 1E-10);

		double biv01 = GaussianUtility.getBivariateCumulativeProbability(x1, x2, true, false, rho);
		assertEquals(0.5864236757641986, biv01, 1E-10);

		double biv10 = GaussianUtility.getBivariateCumulativeProbability(x1, x2, false, true, rho);
		assertEquals(0.10352300242645603, biv10, 1E-10);

		double biv00 = GaussianUtility.getBivariateCumulativeProbability(x1, x2, true, true, rho);
		assertEquals(0.026762827968218002, biv00, 1E-10);
	}
	
	@Test
	public void quantileTests() {
		for (int i = 1; i < 20; i++) {
			double expectedCDFValue = i * .05;
			double quantile = GaussianUtility.getQuantile(expectedCDFValue);
			double cdfValue = GaussianUtility.getCumulativeProbability(quantile);
			assertEquals(expectedCDFValue, cdfValue, 1E-9);
		}
	}
	
	@Test
	public void probabilityDensityTest() {
		double prob = GaussianUtility.getProbabilityDensity(10, 5, 9);
		assertEquals("Comparing probability density", 0.03315905, prob, 1E-8);
		prob = GaussianUtility.getProbabilityDensity(-6, -1, 100);
		assertEquals("Comparing probability density", 0.03520653, prob, 1E-8);
		try {
			prob = GaussianUtility.getProbabilityDensity(-6, -1, -100);
			Assert.fail();	// if get to this line, means that the exception was not throw.
		} catch(InvalidParameterException e) {}
		StandardGaussianDistribution sgd = new StandardGaussianDistribution();
		Matrix mat = new Matrix(1,1);
		mat.setValueAt(0, 0, 1d);
		sgd.getProbabilityDensity(mat);
	}
	
	@Test
	public void testingDensities() {
		double w = 2;
		double mu = 1.2;
		double variance = 2;
		double std = Math.sqrt(variance);
		double dens1 = GaussianUtility.getProbabilityDensity((w - mu) / std) / Math.sqrt(variance);
		double dens2 = GaussianUtility.getProbabilityDensity(w, mu, variance);
		Assert.assertEquals("Comparing density of standardized and non standardized value", dens2, dens1, 1E-8);
	}

}
