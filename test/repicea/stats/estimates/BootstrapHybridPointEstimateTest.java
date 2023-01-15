package repicea.stats.estimates;

import java.util.Random;

import org.junit.Assert;
//import org.junit.Ignore;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.sampling.PopulationUnitWithEqualInclusionProbability;

public class BootstrapHybridPointEstimateTest {

	private static final Random RANDOM = new Random();

	private static int NbRealizations = 25000;
	
	
//	@Ignore
	@Test
	public void simpleTestWithoutModelVariability() {
		PopulationMeanEstimate pe = new PopulationMeanEstimate();
		Matrix obs;
		for (int i = 0; i < 50; i++) {
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, RANDOM.nextGaussian() * 2 + 12);
			pe.addObservation(new PopulationUnitWithEqualInclusionProbability(i + "", obs));
		}

		BootstrapHybridPointEstimate bhpe = new BootstrapHybridPointEstimate(); 
		for (int i = 0; i < 1000; i++) {
			bhpe.addPointEstimate(pe);
		}
		
		double expectedMean = pe.getMean().getValueAt(0, 0);
		double actualMean = bhpe.getMean().getValueAt(0, 0);
		
		System.out.println("Testing without model variability...");
		System.out.println("Expected mean = " + expectedMean + " - actual mean = " + actualMean);
		Assert.assertEquals("Testing mean estimates", expectedMean, actualMean, 1E-8);
		
		
		double expectedVariance = pe.getVariance().getValueAt(0, 0);
		double actualVariance = bhpe.getVariance().getValueAt(0, 0);
		System.out.println("Expected variance = " + expectedVariance + " - actual variance = " + actualVariance);
		
		Assert.assertEquals("Testing variance estimates", expectedVariance, actualVariance, 1E-8);
	}
	

//	@Ignore
	@Test
	public void simpleTestWithoutSamplingVariability() {
		BootstrapHybridPointEstimate bhpe = new BootstrapHybridPointEstimate(); 
		for (int i = 0; i < NbRealizations; i++) {
			PopulationMeanEstimate pe = new PopulationMeanEstimate();
			double deviate = RANDOM.nextGaussian() * 2 + 12;
			Matrix obs;
			for (int j = 0; j < 50; j++) {
				obs = new Matrix(1,1);
				obs.setValueAt(0, 0, deviate);
				pe.addObservation(new PopulationUnitWithEqualInclusionProbability(j + "", obs));
			}
			bhpe.addPointEstimate(pe);
		}
		
		double expectedMean = 12d;
		double actualMean = bhpe.getMean().getValueAt(0, 0);
		System.out.println("Testing without sampling variability...");
		System.out.println("Expected mean = " + expectedMean + " - actual mean = " + actualMean);
		Assert.assertEquals("Testing mean estimates", expectedMean, actualMean, 4E-2);
		
		
		double expectedVariance = 4d;
		double actualVariance = bhpe.getVariance().getValueAt(0, 0);
		System.out.println("Expected variance = " + expectedVariance + " - actual variance = " + actualVariance);
		
		Assert.assertEquals("Testing variance estimates", expectedVariance, actualVariance, 1.5E-1);
	}

//	@Ignore
	@Test
	public void simpleTestWithCompleteVariability() {
		BootstrapHybridPointEstimate bhpe = new BootstrapHybridPointEstimate(); 
		int sampleSize = 10;
		double meanX = 20;
		double stdPopUnit = 15;
		PopulationMeanEstimate pe = new PopulationMeanEstimate();
		Matrix obs;
		for (int j = 0; j < sampleSize; j++) {
			obs = new Matrix(1,1);
			obs.setValueAt(0, 0, meanX + RANDOM.nextGaussian() * stdPopUnit);
			pe.addObservation(new PopulationUnitWithEqualInclusionProbability(j + "", obs));
		}
		
		double mu_x_hat = pe.getMean().getValueAt(0, 0);
		double var_mu_x_hat = pe.getVariance().getValueAt(0, 0);
		
		double meanModel = 0.7;
		double stdModel = .15;
		double stdRes = 1.1;
		
		for (int i = 0; i < NbRealizations; i++) {
			PopulationMeanEstimate peNew = new PopulationMeanEstimate();
			Matrix obsNew;
			double slope = meanModel + RANDOM.nextGaussian() * stdModel; 
			for (String sampleId : pe.getSampleIds()) {
				obsNew = new Matrix(1,1);
				double x = pe.getObservations().get(sampleId).getData().getValueAt(0, 0);
				obsNew.setValueAt(0, 0, x * slope + stdRes * RANDOM.nextGaussian());
				peNew.addObservation(new PopulationUnitWithEqualInclusionProbability(sampleId, obsNew));
			}
			bhpe.addPointEstimate(peNew);
		}
		
		double expectedMean = meanModel * mu_x_hat;
		double actualMean = bhpe.getMean().getValueAt(0, 0);
		
		System.out.println("Testing with complete variability...");
		System.out.println("Expected mean = " + expectedMean + " - actual mean = " + actualMean);
		Assert.assertEquals("Testing mean estimates", expectedMean, actualMean, 0.1);
		
		
		double expectedVariance = mu_x_hat * mu_x_hat * stdModel * stdModel + 
				meanModel * meanModel * var_mu_x_hat -	
				stdModel * stdModel * var_mu_x_hat;	// when dealing with the estimate of the mean, the contribution of the residual error tends to 0, i.e. N * V.bar(e_i) / N^2 = V.bar(e_i) / N. MF2020-12-14
		Matrix varianceMatrix = bhpe.getVariance();
//		System.out.println("Model-related variance = " + varPointEstimate.getModelRelatedVariance());
//		System.out.println("Sampling-related variance = " + varPointEstimate.getSamplingRelatedVariance());
//		System.out.println("Total variance = " + varPointEstimate.getTotalVariance());
		double actualVariance = varianceMatrix.getValueAt(0, 0);

		System.out.println("Expected variance= " + expectedVariance + " - actual variance = " + actualVariance);
		Assert.assertEquals("Testing variance estimates", expectedVariance, actualVariance, 0.8);

		Matrix empiricalCorrection = bhpe.getVarianceBiasCorrection();
		double theoreticalCorrection = -stdModel * stdModel * var_mu_x_hat;
//		System.out.println("Theoretical correction = " + theoreticalCorrection);
//		System.out.println("Empirical correction = " + empiricalCorrection);
		Assert.assertEquals("Comparing variance bias correction", theoreticalCorrection, empiricalCorrection.getValueAt(0, 0), 2E-2);
	}

	@Test
	public void multivariateTestWithCompleteVariability() {
		BootstrapHybridPointEstimate bhpe = new BootstrapHybridPointEstimate(); 
		int sampleSize = 10;
		Matrix meanX = new Matrix(2,1);
		meanX.setValueAt(0, 0, 20d);
		meanX.setValueAt(1, 0, 12d);
		SymmetricMatrix varianceX = new SymmetricMatrix(2);
		varianceX.setValueAt(0, 0, 10d * 10d);
		varianceX.setValueAt(1, 1, 4d * 4d);
		varianceX.setValueAt(1, 0, 4d * 10d * .5);

		GaussianEstimate popUnitGenerator = new GaussianEstimate(meanX, varianceX);
		PopulationMeanEstimate pe = new PopulationMeanEstimate();
		for (int j = 0; j < sampleSize; j++) {
			Matrix randomObs = popUnitGenerator.getRandomDeviate();
			pe.addObservation(new PopulationUnitWithEqualInclusionProbability(j + "", randomObs));
		}
		
		Matrix mu_x_hat = pe.getMean();
		Matrix var_mu_x_hat = pe.getVariance();
		
		double meanModel = 0.7;
		double stdModel = .15;
		double stdRes = 1.1;
		
		for (int i = 0; i < NbRealizations; i++) {
			PopulationMeanEstimate peNew = new PopulationMeanEstimate();
			Matrix obsNew;
			Matrix slope = new Matrix(1,1,meanModel + RANDOM.nextGaussian() * stdModel,0); 
			for (String sampleId : pe.getSampleIds()) {
				Matrix x = pe.getObservations().get(sampleId).getData();
				obsNew =  x.multiply(slope).add(StatisticalUtility.drawRandomVector(x.m_iRows, Distribution.Type.GAUSSIAN).scalarMultiply(stdRes));
				peNew.addObservation(new PopulationUnitWithEqualInclusionProbability(sampleId, obsNew));
			}
			bhpe.addPointEstimate(peNew);
		}
		
		Matrix expectedMean = mu_x_hat.scalarMultiply(meanModel);
		Matrix actualMean = bhpe.getMean();


		System.out.println("Testing multivariate with complete variability...");
		System.out.println("Expected mean = " + expectedMean + " - actual mean = " + actualMean);
		Assert.assertTrue("Testing mean estimates", !expectedMean.subtract(actualMean).getAbsoluteValue().anyElementLargerThan(0.1));
		
		
		Matrix expectedVariance = mu_x_hat.multiply(mu_x_hat.transpose()).scalarMultiply(stdModel * stdModel)
				.add(var_mu_x_hat.scalarMultiply(meanModel * meanModel))
				.subtract(var_mu_x_hat.scalarMultiply(stdModel * stdModel));	// when dealing with the estimate of the mean, the contribution of the residual error tends to 0, i.e. N * V.bar(e_i) / N^2 = V.bar(e_i) / N. MF2020-12-14

//		System.out.println("Model-related variance = " + varPointEstimate.getModelRelatedVariance());
//		System.out.println("Sampling-related variance = " + varPointEstimate.getSamplingRelatedVariance());
		Matrix actualVariance = bhpe.getVariance();
//		System.out.println("Total variance = " + actualVariance);
//		System.out.println("Expected variance = " + expectedVariance);

		System.out.println("Expected variance = " + expectedVariance + " - actual variance = " + actualVariance);
		Assert.assertTrue("Testing variance estimates", !expectedVariance.subtract(actualVariance).getAbsoluteValue().anyElementLargerThan(3E-1));
		
		Matrix empiricalCorrection = bhpe.getVarianceBiasCorrection();
		Matrix theoreticalCorrection = var_mu_x_hat.scalarMultiply(-stdModel * stdModel);
//		System.out.println("Theoretical correction = " + theoreticalCorrection);
//		System.out.println("Empirical correction = " + empiricalCorrection);
		Matrix relDiff = empiricalCorrection.elementWiseDivide(theoreticalCorrection).scalarAdd(-1d);

		Assert.assertTrue("Testing variance estimates", !relDiff.getAbsoluteValue().anyElementLargerThan(2E-1));
	}
	
	
}
