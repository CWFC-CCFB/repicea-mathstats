/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Wood Fibre Centre, Canadian Forest Service
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 3 of the License, or (at your option) any later version.
 *
 * This library is distributed with the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied
 * warranty of MERCHANTABILITY or FITNESS FOR A
 * PARTICULAR PURPOSE. See the GNU Lesser General Public
 * License for more details.
 *
 * Please see the license at http://www.gnu.org/copyleft/lesser.html.
 */
package repicea.math;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.lang.REpiceaSystem;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.ComplexMonteCarloEstimate;
import repicea.stats.estimates.MonteCarloEstimate;

/**
 * The ComplexNumberSimpleCaseStudy class implements a simulation study around 
 * the exponential, square and log transformation of Gaussian variables.
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberArsinhTransformationBalancedSimulationStudy extends AbstractComplexNumberSimulationStudy {


	/**
	 * A model fitted to a sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	static class Model extends AbstractComplexNumberSimulationStudy.Model {
		
		Model(Sample sample, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			super(sample, betaHat, invXtX, sigma2Hat, upsilon);
		}

		private Matrix getRandomDeviate(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			double sigma2_b = sigma2Hat + sigma2Hat * (chiSquareDeviate/upsilon - 1);
			Matrix betaHat_b = betaHat.add(
					getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(Math.sqrt(sigma2_b)));
			Matrix xBeta_b = xMat.multiply(betaHat_b);
			Matrix sinhResult = xBeta_b.expMatrix().subtract(xBeta_b.scalarMultiply(-1d).expMatrix()).scalarMultiply(0.5);
			return sinhResult.scalarMultiply(Math.exp(sigma2_b * .5));
		}
		
		private Matrix getNaiveEstimator(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues);
			Matrix sinhResult = xBeta.expMatrix().subtract(xBeta.scalarMultiply(-1d).expMatrix()).scalarMultiply(0.5);
			return sinhResult.scalarMultiply(Math.exp(sigma2Hat * .5));
		}
		
		private ComplexMatrix getComplexRandomDeviate(List<Double> xValues) {
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate/upsilon - 1));
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();
			Matrix xMat = createMatrixX(xValues);
			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)); // here we get C * epsilon
			
			ComplexNumber[] cnArrayMean = new ComplexNumber[xValues.size()];
			for (int ii = 0; ii < xValues.size(); ii++) {
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0)); // here we get 0 + x C epsilon i 
				betaDeviates = betaDeviates.multiply(sigmaHat_b); // here we get 0 + sqrt{sigma2Hat_b} x C epsilon i
				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
				ComplexNumber xBeta_b = betaDeviates.add(xBetaHat); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i
				ComplexNumber sinhResult = xBeta_b.exp().subtract(xBeta_b.multiply(-1d).exp()).multiply(.5);
				ComplexNumber sinhResultPlusExpSigma2DividedBy2 = sinhResult.multiply(sigma2Hat_b.multiply(.5).exp()); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i + 0.5 * sigma2_b  				
				cnArrayMean[ii] = sinhResultPlusExpSigma2DividedBy2;
			}
			return new ComplexMatrix(cnArrayMean); 
		}
	}
	
	private void doRun(int sampleSize, 
			int nbRealizations, 
			int nbCMCRealization,
			double b0, 
			double b1, 
			double s2,
			String filename) throws IOException {
		System.out.println("Simulating [" + b0 + "; " + b1 + "; " + s2 +"] with sample size n = " + sampleSize);
		CSVWriter writer = new CSVWriter(new File(filename.concat("" + sampleSize + ".csv")), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("x"));
		fields.add(new CSVField("Naive"));
		fields.add(new CSVField("MC_meanEst"));
		fields.add(new CSVField("CMC_meanEst_real"));
		fields.add(new CSVField("CMC_meanEst_imag"));
		fields.add(new CSVField("CMC_varEst_real"));
		fields.add(new CSVField("CMC_varEst_imag"));
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		trueBeta.setValueAt(0, 0, b0);
		trueBeta.setValueAt(1, 0, b1);
		double trueVariance = s2;
		
		
		List<Double> xValues = Sample.X_VALUES;
		for (int real = 1; real <= nbRealizations; real++) {
			if (real % 1000 == 0) {
				System.out.println("Running realization " + real);
			}

			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = (Model) s.getModel(this);
			
			Matrix naiveEstimator = m.getNaiveEstimator(xValues);
			
			MonteCarloEstimate mcEstimator = new MonteCarloEstimate();
			ComplexMonteCarloEstimate cmcEstimator = new ComplexMonteCarloEstimate();

			for (int innerReal = 0; innerReal < nbCMCRealization; innerReal++) {
				mcEstimator.addRealization(m.getRandomDeviate(xValues)); 
				ComplexMatrix complexRealizations = m.getComplexRandomDeviate(xValues);
				cmcEstimator.addRealization(complexRealizations); 
			}
			
			Matrix mcMean = mcEstimator.getMean();
			
			ComplexMatrix cmcMean = cmcEstimator.getMean();
			ComplexSymmetricMatrix cmcPsVar = cmcEstimator.getPseudoVariance();

			for (int ii = 0; ii < xValues.size(); ii++) {
				Object[] record = new Object[fields.size()];
				record[0] = real;
				record[1] = m.betaHat.getValueAt(0, 0);
				record[2] = m.betaHat.getValueAt(1, 0);
				record[3] = m.sigma2Hat;
				record[4] = xValues.get(ii);
				record[5] = naiveEstimator.getValueAt(ii, 0);
				record[6] = mcMean.getValueAt(ii, 0);
				record[7] = cmcMean.getValueAt(ii, 0).realPart;
				record[8] = cmcMean.getValueAt(ii, 0).imaginaryPart;
				record[9] = cmcPsVar.getValueAt(ii, ii).realPart;
				record[10] = cmcPsVar.getValueAt(ii, ii).imaginaryPart;
				writer.addRecord(record);
			}
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException, InterruptedException {
		String filename = REpiceaSystem.retrieveArgument("-outdir", Arrays.asList(args));
		System.out.println("Export filename: " + filename);
		List<Double> variances = new ArrayList<Double>();
		variances.add(1.0);
		variances.add(2.0);
		for (double variance : variances) {
			String rootFilename = filename.concat("caseStudySinhBalanced_" + variance + "_");
			List<Integer> sampleSizes = new ArrayList<Integer>();
			sampleSizes.add(25);
			sampleSizes.add(50);
			sampleSizes.add(100);
			sampleSizes.add(200);
			List<Thread> threads = new ArrayList<Thread>();
			for (int sampleSize : sampleSizes) {
				Runnable doRun = new Runnable() {
					public void run() {
						ComplexNumberArsinhTransformationBalancedSimulationStudy s = new ComplexNumberArsinhTransformationBalancedSimulationStudy();
						try {
							s.doRun(sampleSize, 50000, 10000, 1, 0.25, variance, rootFilename);
						} catch (IOException e) {
							e.printStackTrace();
						}
					};
				};
				Thread t = new Thread(doRun);
				threads.add(t);
				t.start();
			}
			for (Thread t : threads) {
				t.join();
			}
		}
	}

	@Override
	ComplexNumberArsinhTransformationBalancedSimulationStudy.Model createModel(Sample s, 
			Matrix betaHat, 
			SymmetricMatrix invXtX, 
			double sigma2Hat, 
			int upsilon) {
		return new ComplexNumberArsinhTransformationBalancedSimulationStudy.Model(s, 
				betaHat, 
				invXtX, 
				sigma2Hat, 
				upsilon);
	}	
}
