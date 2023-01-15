package repicea.stats.estimates;

import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.sampling.PopulationUnitWithUnequalInclusionProbability;

public class HorvitzThompsonTauEstimateTest {

	@Test
	public void simpleTotalAndVarianceEstimateTest() {
		int populationSize = 1000;
		List<Double> sample = new ArrayList<Double>();
		sample.add(2d);
		sample.add(4d);
		sample.add(2d);
		sample.add(5d);
		sample.add(7d);
		sample.add(1d);
		sample.add(5d);
		sample.add(4d);
		sample.add(7d);
		PopulationTotalEstimate estimate = new PopulationTotalEstimate();
		Matrix obs;
		for (int i = 0; i < sample.size(); i++) {
			double value = sample.get(i);
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, value);
			estimate.addObservation(new PopulationUnitWithUnequalInclusionProbability(i + "", obs, 1d/populationSize));
		}
		Matrix total = estimate.getMean();
		Assert.assertEquals("Testing the estimate of the total", 4111.11111111111, total.getValueAt(0, 0), 1E-8);
		Matrix totalVariance = estimate.getVariance();
		Assert.assertEquals("Testing the variance of the total", 507734.5679, totalVariance.getValueAt(0, 0), 1E-4);
	}
	
	
	
	
}
