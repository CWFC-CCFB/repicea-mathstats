package repicea.stats.distributions;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;

public class UniformDistributionTest {

	@Test
	public void simpleMeanTest() {
		Matrix lowerBound = new Matrix(2,1);
		lowerBound.setValueAt(0, 0, -2d);
		lowerBound.setValueAt(1, 0, 1d);
			
		Matrix upperBound = new Matrix(2,1);
		upperBound.setValueAt(0, 0, 1d);
		upperBound.setValueAt(1, 0, 3d);
		
		UniformDistribution dist = new UniformDistribution(lowerBound, upperBound);
		Matrix mean = dist.getMean();
		Assert.assertEquals(-.5, mean.getValueAt(0, 0), 1E-8);
		Assert.assertEquals(2d, mean.getValueAt(1, 0), 1E-8);
	}
	
	
	@Test
	public void simpleVarianceTest() {
		Matrix lowerBound = new Matrix(2,1);
		lowerBound.setValueAt(0, 0, -2d);
		lowerBound.setValueAt(1, 0, 1d);
			
		Matrix upperBound = new Matrix(2,1);
		upperBound.setValueAt(0, 0, 1d);
		upperBound.setValueAt(1, 0, 3d);
		
		UniformDistribution dist = new UniformDistribution(lowerBound, upperBound);
		Matrix variance = dist.getVariance();
		Assert.assertEquals(.75, variance.getValueAt(0, 0), 1E-8);
		Assert.assertEquals(1d/3, variance.getValueAt(1, 1), 1E-8);
	}

	@Test
	public void simpleMonteCarloRandomDeviateTest() {
		Matrix lowerBound = new Matrix(2,1);
		lowerBound.setValueAt(0, 0, -2d);
		lowerBound.setValueAt(1, 0, 1d);
			
		Matrix upperBound = new Matrix(2,1);
		upperBound.setValueAt(0, 0, 1d);
		upperBound.setValueAt(1, 0, 3d);
		
		UniformDistribution dist = new UniformDistribution(lowerBound, upperBound);
		int nbMc = 1000000;
		Matrix meanMC = new Matrix(2,1); 
		for (int i = 0; i < nbMc; i++) {
			meanMC = meanMC.add(dist.getRandomRealization());
		}
		meanMC = meanMC.scalarMultiply(1d/nbMc);
		Matrix mean = dist.getMean();
		Matrix relativeResult = meanMC.elementWiseDivide(mean);
		Assert.assertEquals(1d, relativeResult.getValueAt(0, 0), 5E-3);
		Assert.assertEquals(1d, relativeResult.getValueAt(1, 0), 5E-3);
	}


	@Test
	public void simpleDensityTest() {
		Matrix lowerBound = new Matrix(2,1);
		lowerBound.setValueAt(0, 0, -2d);
		lowerBound.setValueAt(1, 0, 1d);
			
		Matrix upperBound = new Matrix(2,1);
		upperBound.setValueAt(0, 0, 1d);
		upperBound.setValueAt(1, 0, 3d);
		
		UniformDistribution dist = new UniformDistribution(lowerBound, upperBound);
		Matrix values = new Matrix(2,1,0,2);
		double actualDensity = dist.getProbabilityDensity(values);
		double expectedDensity = 1d/(upperBound.getValueAt(0, 0) - lowerBound.getValueAt(0, 0)) * 1d/(upperBound.getValueAt(1, 0) - lowerBound.getValueAt(1, 0));
		Assert.assertEquals("Testing density", expectedDensity, actualDensity, 1E-16);
	}

	@Test
	public void simpleDensityTestBeyondBounds() {
		Matrix lowerBound = new Matrix(2,1);
		lowerBound.setValueAt(0, 0, -2d);
		lowerBound.setValueAt(1, 0, 1d);
			
		Matrix upperBound = new Matrix(2,1);
		upperBound.setValueAt(0, 0, 1d);
		upperBound.setValueAt(1, 0, 3d);
		
		UniformDistribution dist = new UniformDistribution(lowerBound, upperBound);
		Matrix values = new Matrix(2,1,0,4);
		double actualDensity = dist.getProbabilityDensity(values);
		double expectedDensity = 0d;
		Assert.assertEquals("Testing density", expectedDensity, actualDensity, 1E-16);
	}

	
}
