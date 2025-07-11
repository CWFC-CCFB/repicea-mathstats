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
public class ComplexNumberLogTransformationBalancedSimulationStudy extends AbstractComplexNumberSimulationStudy {


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
			Matrix result = xMat.multiply(betaHat_b).scalarAdd(0.5 * sigma2_b);
			return result.expMatrix();
		}
		
		private Matrix getBeauchampAndOlsonEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix xBeta = getXBeta(xValues);
			
			Matrix output = new Matrix(xValues.size(), 1);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				double term1 = Math.exp(xBeta.getValueAt(ii, 0) + sigma2Hat * .5);
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				int n = upsilon + xMat_ii.m_iCols;
				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
				double phi = xInvXtXxT.getValueAt(0, 0) * n;
				double sigma2hat2 = sigma2Hat * sigma2Hat;
				double term2 = 1 - sigma2Hat * (2*phi + sigma2Hat)/(4*n) + 
						sigma2hat2 * (sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2Hat + 4 * phi * phi + 16 * phi) / (32 * n * n);
				output.setValueAt(ii, 0, term1 * term2);
			}
			return output;
		}

		private Matrix getBeauchampAndOlsonVarianceEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix xBeta = getXBeta(xValues);

			Matrix expPart = xBeta.scalarMultiply(2d).scalarAdd(sigma2Hat).expMatrix();
			
			Matrix output = new Matrix(xValues.size(), 1);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				int n = upsilon + xMat_ii.m_iCols;
				Matrix varXBetaMat = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose()).scalarMultiply(sigma2Hat);
				double varXBeta = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose()).scalarMultiply(sigma2Hat).getValueAt(0, 0);
				double sigma2hat2Term = sigma2Hat * sigma2Hat / (2*n);
				output.setValueAt(ii, 0, expPart.getValueAt(ii, 0) * (varXBeta + sigma2hat2Term));
			}
			return output;
		}

		private Matrix getBaskervilleEstimator(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues);
			return xBeta.scalarAdd(sigma2Hat * .5).expMatrix();
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
				ComplexNumber mean = betaDeviates.add(xBetaHat); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i
				ComplexNumber meanPlusCorrectionFactor = mean.add(sigma2Hat_b.multiply(0.5)); // here we get x beta + + sqrt{sigma2Hat_b} x C epsilon i + 0.5 * sigma2_b  				
				cnArrayMean[ii] = meanPlusCorrectionFactor.exp();
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
		fields.add(new CSVField("Baskerville"));
		fields.add(new CSVField("Beauchamp"));
		fields.add(new CSVField("BeauchampVariance"));
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
			
			Matrix beauchampAndOlsonEstimator = m.getBeauchampAndOlsonEstimator(xValues);
			Matrix beauchampAndOlsonVarianceEstimator = m.getBeauchampAndOlsonVarianceEstimator(xValues);
			Matrix baskervilleEstimator = m.getBaskervilleEstimator(xValues);
			
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
				record[5] = baskervilleEstimator.getValueAt(ii, 0);
				record[6] = beauchampAndOlsonEstimator.getValueAt(ii, 0);
				record[7] = beauchampAndOlsonVarianceEstimator.getValueAt(ii, 0);
				record[8] = mcMean.getValueAt(ii, 0);
				record[9] = cmcMean.getValueAt(ii, 0).realPart;
				record[10] = cmcMean.getValueAt(ii, 0).imaginaryPart;
				record[11] = cmcPsVar.getValueAt(ii, ii).realPart;
				record[12] = cmcPsVar.getValueAt(ii, ii).imaginaryPart;
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
			String rootFilename = filename.concat("caseStudyBalanced_" + variance + "_");
			List<Integer> sampleSizes = new ArrayList<Integer>();
			sampleSizes.add(25);
			sampleSizes.add(50);
			sampleSizes.add(100);
			sampleSizes.add(200);
			List<Thread> threads = new ArrayList<Thread>();
			for (int sampleSize : sampleSizes) {
				Runnable doRun = new Runnable() {
					public void run() {
						ComplexNumberLogTransformationBalancedSimulationStudy s = new ComplexNumberLogTransformationBalancedSimulationStudy();
						try {
							s.doRun(sampleSize, 50000, 10000, 2, 0.25, variance, rootFilename);
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
	ComplexNumberLogTransformationBalancedSimulationStudy.Model createModel(Sample s, 
			Matrix betaHat, 
			SymmetricMatrix invXtX, 
			double sigma2Hat, 
			int upsilon) {
		return new ComplexNumberLogTransformationBalancedSimulationStudy.Model(s, 
				betaHat, 
				invXtX, 
				sigma2Hat, 
				upsilon);
	}	
}
