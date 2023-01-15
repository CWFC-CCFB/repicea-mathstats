package repicea.stats.estimates;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.math.Matrix;
import repicea.math.SymmetricMatrix;
import repicea.util.ObjectUtility;

public class ProductOfEstimates {

	private static double LowFactor = 0.05;
	private static double HighFactor = 0.20;
	private static enum Method {Naive, Propagation}


	private static SimpleEstimate getEstimate(List<Estimate> estimates, Method method) {
		Estimate currentEstimate = null;
		for (int i = 1; i < estimates.size(); i++) {
			if (i == 1) {
				currentEstimate = estimates.get(i - 1);
			}

			Matrix alphaMean = currentEstimate.getMean();
			Matrix betaMean = estimates.get(i).getMean();
			SymmetricMatrix alphaVariance = currentEstimate.getVariance();
			SymmetricMatrix betaVariance = estimates.get(i).getVariance();
			Matrix newMean = alphaMean.multiply(betaMean);
			SymmetricMatrix newVariance = SymmetricMatrix.convertToSymmetricIfPossible( 
					alphaMean.elementWisePower(2d).multiply(betaVariance).add(betaMean.elementWisePower(2d).multiply(alphaVariance))
					);
			if (method == Method.Naive) {
				newVariance = (SymmetricMatrix) newVariance.add(alphaVariance.multiply(betaVariance));
			}
			if (newVariance.getValueAt(0, 0) < 0d) {
				throw new UnsupportedOperationException("The product of the estimates yields negative variance");
			}
			currentEstimate = new SimpleEstimate(newMean, newVariance);
		}
		return (SimpleEstimate) currentEstimate;
	}

	private static MonteCarloEstimate getMonteCarloEstimate(List<Estimate> estimates, int nbRealizations) {
		MonteCarloEstimate est = new MonteCarloEstimate(); 
		for (int real = 0; real < nbRealizations; real++) {
			double currentProduct = -1;
			for (int i = 1; i < estimates.size(); i++) {
				if (i == 1) {
					currentProduct = estimates.get(i - 1).getRandomDeviate().getValueAt(0, 0);
				}

				currentProduct *= estimates.get(i).getRandomDeviate().getValueAt(0, 0);
			}
			Matrix obs = new Matrix(1,1);
			obs.setValueAt(0, 0, currentProduct);
			est.addRealization(obs);
		}
		return est;
	}

	private static Estimate getEstimate(Estimate trueMean, VarianceEstimate trueVariance, boolean isLogNormal) {
		double mean = trueMean.getRandomDeviate().getValueAt(0, 0);
		double variance = trueVariance.getRandomDeviate().getValueAt(0, 0);
		if (isLogNormal) {
			return new LogNormalEstimate(mean, variance, false);
		} else {
			return new GaussianEstimate(mean, variance);
		}
	}

	private static VarianceEstimate getTrueVariance(int df, double mean, boolean lowVariability) {
		if (lowVariability) {
			return new VarianceEstimate(df, LowFactor * LowFactor * mean * mean);
		} else {
			return new VarianceEstimate(df, HighFactor * HighFactor * mean * mean);
		}
	}


	private static void runSimulation(int nbMaxRealization,
									  boolean lowAlpha,
									  boolean lowBeta,
									  boolean lowGamma,
									  boolean useLogNormal) throws Exception {
		runSimulation(nbMaxRealization,
				lowAlpha,
				lowBeta,
				lowGamma,
				0d,
				0d,
				0d,
				useLogNormal);
	}

	private static void runSimulation(int nbMaxRealization,
									  boolean lowAlpha,
									  boolean lowBeta,
									  boolean lowGamma,
									  double biasAlpha,
									  double biasBeta,
									  double biasGamma,
									  boolean useLogNormal) throws Exception {
		String simulationName = getFilenameSuffix(lowAlpha, lowBeta, lowGamma, biasAlpha == 0d && biasBeta == 0d && biasGamma == 0d, useLogNormal);
		System.out.println("Running simulation with lowAlpha = " + ((Boolean) lowAlpha).toString() + " lowBeta = " + ((Boolean) lowBeta).toString() + " lowGamma = " + ((Boolean) lowGamma).toString() + " ...");
		int df = 300;
		double alpha = 20d;
		double beta = 10d;
		double gamma = 2d;
		double trueMean = alpha * beta * gamma;

		VarianceEstimate trueVarAlpha = ProductOfEstimates.getTrueVariance(df, alpha, lowAlpha);
		VarianceEstimate trueVarBeta = ProductOfEstimates.getTrueVariance(df, beta, lowBeta);
		VarianceEstimate trueVarGamma = ProductOfEstimates.getTrueVariance(df, gamma, lowGamma);
		double varAlpha = trueVarAlpha.getMean().getValueAt(0, 0);
		double varBeta = trueVarBeta.getMean().getValueAt(0, 0);
		double varGamma = trueVarGamma.getMean().getValueAt(0, 0);

		double biasedAlpha = alpha * (1d + biasAlpha);
		double biasedBeta = beta * (1d + biasBeta);
		double biasedGamma = gamma * (1d + biasGamma);

		Estimate expectedAlpha;
		Estimate expectedBeta;
		Estimate expectedGamma;
		if (useLogNormal) {
			expectedAlpha = new LogNormalEstimate(biasedAlpha, varAlpha, false);
			expectedBeta = new LogNormalEstimate(biasedBeta, varBeta, false);
			expectedGamma = new LogNormalEstimate(biasedGamma, varGamma, false);
		} else {
			expectedAlpha = new GaussianEstimate(biasedAlpha, varAlpha);
			expectedBeta = new GaussianEstimate(biasedBeta, varBeta);
			expectedGamma = new GaussianEstimate(biasedGamma, varGamma);
		}

		MonteCarloEstimate muGoodman = new MonteCarloEstimate();		
		MonteCarloEstimate varGoodman = new MonteCarloEstimate();		
		MonteCarloEstimate muNaive = new MonteCarloEstimate();		
		MonteCarloEstimate varNaive = new MonteCarloEstimate();		
		MonteCarloEstimate muPropagation = new MonteCarloEstimate();		
		MonteCarloEstimate varPropagation = new MonteCarloEstimate();		
		MonteCarloEstimate muMonteCarlo = new MonteCarloEstimate();		
		MonteCarloEstimate varMonteCarlo = new MonteCarloEstimate();		
		MonteCarloEstimate coverage = new MonteCarloEstimate();
		MonteCarloEstimate muRescaledMonteCarlo = new MonteCarloEstimate();		
		MonteCarloEstimate varRescaledMonteCarlo = new MonteCarloEstimate();		
		MonteCarloEstimate rescaledCoverage = new MonteCarloEstimate();
		MonteCarloEstimate mse = new MonteCarloEstimate();

		for (int real = 0; real < nbMaxRealization; real++) {
			if (real%1000 == 0) {
				System.out.println("Realization " + real);
			}
			List<Estimate> estimates = new ArrayList<Estimate>();
			estimates.add(getEstimate(expectedAlpha, trueVarAlpha, useLogNormal));
			estimates.add(getEstimate(expectedBeta, trueVarBeta, useLogNormal));
			estimates.add(getEstimate(expectedGamma, trueVarGamma, useLogNormal));

			SimpleEstimate productGoodman = Estimate.getProductOfManyEstimates(estimates);
			muGoodman.addRealization(productGoodman.getMean());
			varGoodman.addRealization(productGoodman.getVariance());

			SimpleEstimate productNaive = getEstimate(estimates, Method.Naive);
			muNaive.addRealization(productNaive.getMean());
			varNaive.addRealization(productNaive.getVariance());

			SimpleEstimate productPropagation = getEstimate(estimates, Method.Propagation);
			muPropagation.addRealization(productPropagation.getMean());
			varPropagation.addRealization(productPropagation.getVariance());

			MonteCarloEstimate productMC = getMonteCarloEstimate(estimates, 5000);
			muMonteCarlo.addRealization(productMC.getMean());
			varMonteCarlo.addRealization(productMC.getVariance());
			ConfidenceInterval ci = productMC.getConfidenceIntervalBounds(.95);
			Matrix in = new Matrix(1,1);
			if (ci.getLowerLimit().getValueAt(0, 0) <= trueMean && trueMean <= ci.getUpperLimit().getValueAt(0, 0)) {
				in.setValueAt(0, 0, 1d);
			}
			coverage.addRealization(in);


			MonteCarloEstimate scaledProductMC = new MonteCarloEstimate();
			double scalingFactor = Math.sqrt(productGoodman.getVariance().getValueAt(0, 0) / productNaive.getVariance().getValueAt(0, 0));
			if (Double.isNaN(scalingFactor)) {
				throw new UnsupportedOperationException("Trying to compute a negative square root when rescaling the Monte Carlo variance estimator");
			}
			double mean = productMC.getMean().getValueAt(0, 0);
			Matrix newReal;
			for (Matrix r : productMC.getRealizations()) {
				newReal = new Matrix(1,1);
				double rValue = r.getValueAt(0, 0);
				newReal.setValueAt(0, 0, scalingFactor * (rValue - mean) + mean);
				scaledProductMC.addRealization(newReal);
			}
			muRescaledMonteCarlo.addRealization(scaledProductMC.getMean());
			varRescaledMonteCarlo.addRealization(scaledProductMC.getVariance());


			ci = scaledProductMC.getConfidenceIntervalBounds(.95);
			in = new Matrix(1,1);
			if (ci.getLowerLimit().getValueAt(0, 0) <= trueMean && trueMean <= ci.getUpperLimit().getValueAt(0, 0)) {
				in.setValueAt(0, 0, 1d);
			}
			rescaledCoverage.addRealization(in);

			Matrix error = productGoodman.getMean().scalarAdd(- trueMean);
			Matrix se = error.transpose().multiply(error);
			mse.addRealization(se);

		}

		String filename = ObjectUtility.getPackagePath(ProductOfEstimates.class).replace("bin", "test").concat("simulation" + simulationName + ".csv");
		CSVWriter writer = new CSVWriter(new File(filename), false, ","); // splitter is now "," and not ";"
		List<FormatField> formatFields = new ArrayList<FormatField>();
		formatFields.add(new CSVField("sampleSize"));
		formatFields.add(new CSVField("trueMean"));
		formatFields.add(new CSVField("empMean"));
		formatFields.add(new CSVField("trueVariance"));
		formatFields.add(new CSVField("empVariance"));
		formatFields.add(new CSVField("goodmanVariance"));
		formatFields.add(new CSVField("naiveVariance"));
		formatFields.add(new CSVField("propagationVariance"));
		formatFields.add(new CSVField("MeanMC"));
		formatFields.add(new CSVField("VarMC"));
		formatFields.add(new CSVField("coverage"));
		formatFields.add(new CSVField("MeanRescaledMC"));
		formatFields.add(new CSVField("VarRescaledMC"));
		formatFields.add(new CSVField("rescaledCoverage"));
		formatFields.add(new CSVField("mse"));

		writer.setFields(formatFields);


		Object[] record = new Object[15];
		record[0] = simulationName;
		record[1] = trueMean;
		record[2] = muGoodman.getMean().getValueAt(0, 0);
		double expMeanAB = biasedAlpha * biasedBeta;
		double trueVarianceAB = biasedAlpha * biasedAlpha * varBeta + varAlpha * biasedBeta * biasedBeta + varAlpha * varBeta; 
		double trueVariance = expMeanAB * expMeanAB * varGamma + trueVarianceAB * biasedGamma * biasedGamma + trueVarianceAB * varGamma; 
		record[3] = trueVariance;
		record[4] = muGoodman.getVariance().getValueAt(0, 0);
		record[5] = varGoodman.getMean().getValueAt(0, 0);
		record[6] = varNaive.getMean().getValueAt(0, 0);
		record[7] = varPropagation.getMean().getValueAt(0, 0);
		record[8] = muMonteCarlo.getMean().getValueAt(0, 0);
		record[9] = varMonteCarlo.getMean().getValueAt(0, 0);
		record[10] = coverage.getMean().getValueAt(0, 0);
		record[11] = muRescaledMonteCarlo.getMean().getValueAt(0, 0);
		record[12] = varRescaledMonteCarlo.getMean().getValueAt(0, 0);
		record[13] = rescaledCoverage.getMean().getValueAt(0, 0);
		record[14] = mse.getMean().getValueAt(0, 0);

		writer.addRecord(record);
		writer.close();

	}

	private static String getFilenameSuffix(boolean b1, boolean b2, boolean b3, boolean unbiased, boolean useLogNormals) {
		String suffix = "";
		if (b1) {
			suffix += "H";
		} else {
			suffix += "L";
		}
		if (b2) {
			suffix += "H";
		} else {
			suffix += "L";
		}
		if (b3) {
			suffix += "H";
		} else {
			suffix += "L";
		}
		if (!unbiased) {
			suffix += "_b";
		}
		if (useLogNormals) {
			suffix += "_ln";
		}
		return suffix;
	}
	
	public static void main(String[] args) throws Exception {
		int nbRealizations = 50000;

		// Using Gaussian distributions
//		runSimulation(nbRealizations, true, true, true, false);
//		runSimulation(nbRealizations, false, true, true, false);
//		runSimulation(nbRealizations, true, false, true, false);
//		runSimulation(nbRealizations, true, true, false, false);
//		runSimulation(nbRealizations, false, false, true, false);
//		runSimulation(nbRealizations, false, true, false, false);
//		runSimulation(nbRealizations, true, false, false, false);
//		runSimulation(nbRealizations, false, false, false, false);


//		runSimulation(nbRealizations, true, true, true, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, false, true, true, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, true, false, true, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, true, true, false, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, false, false, true, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, false, true, false, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, true, false, false, 0.02, 0.02, 0.02, false);
//		runSimulation(nbRealizations, false, false, false, 0.02, 0.02, 0.02, false);

		// Using log normal distributions
		runSimulation(nbRealizations, true, true, true, true);
		runSimulation(nbRealizations, false, true, true, true);
		runSimulation(nbRealizations, true, false, true, true);
		runSimulation(nbRealizations, true, true, false, true);
		runSimulation(nbRealizations, false, false, true, true);
		runSimulation(nbRealizations, false, true, false, true);
		runSimulation(nbRealizations, true, false, false, true);
		runSimulation(nbRealizations, false, false, false, true);


		runSimulation(nbRealizations, true, true, true, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, false, true, true, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, true, false, true, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, true, true, false, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, false, false, true, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, false, true, false, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, true, false, false, 0.02, 0.02, 0.02, true);
		runSimulation(nbRealizations, false, false, false, 0.02, 0.02, 0.02, true);


	}

}
