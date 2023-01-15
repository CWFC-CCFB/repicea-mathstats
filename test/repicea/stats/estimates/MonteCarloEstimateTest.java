package repicea.stats.estimates;

import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.sampling.PopulationUnitWithEqualInclusionProbability;

public class MonteCarloEstimateTest {

	private static int NbRealizations = 50000;
	
	@Test
	public void percentilesTest() {
		MonteCarloEstimate estimate = new MonteCarloEstimate();
		Random random = new Random();
		Matrix m;
		for (int i = 0; i < NbRealizations; i++) {
			m = new Matrix(2,1);
			m.setValueAt(0, 0, random.nextDouble());
			m.setValueAt(1, 0, random.nextDouble() * 2);
			estimate.addRealization(m);
		}
		
		Matrix percentiles = estimate.getQuantileForProbability(0.95);
		Assert.assertEquals(0.95, percentiles.getValueAt(0, 0), 1E-2);
		Assert.assertEquals(0.95 * 2, percentiles.getValueAt(1, 0), 1E-2);
	}
	
	@Test
	public void univariateTotalVarianceLawTest() {
		double outerStd = 5d;
		double innerStd = 5d;
		int sampleSize = 100;
		Random generator = new Random();
		LawOfTotalVarianceMonteCarloEstimate output = new LawOfTotalVarianceMonteCarloEstimate();
		for (int i = 0; i < 10000; i++) {
			double meanForThisRealization = generator.nextGaussian() * outerStd;
			PopulationMeanEstimate estimate = new PopulationMeanEstimate();
			Matrix obs;
			for (int j = 0; j < sampleSize; j++) {
				obs = new Matrix(1,1);
				obs.setValueAt(0, 0, generator.nextGaussian() * innerStd + meanForThisRealization); 
				estimate.addObservation(new PopulationUnitWithEqualInclusionProbability(j + "", obs));
			}
			output.addRealization(estimate);
		}
		double totalVariance = output.getVariance().getValueAt(0, 0);
		double expectedVariance = outerStd * outerStd + innerStd * innerStd / sampleSize;
		double relativeDifference = (totalVariance - expectedVariance)/expectedVariance; 
		System.out.println("Relative difference of " + relativeDifference);
		Assert.assertEquals(0d, relativeDifference, 0.4);
	}
	
}
