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
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.stats.Distribution;
import repicea.stats.Distribution.Type;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.ComplexMonteCarloEstimate;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.util.ObjectUtility;

/**
 * The ComplexNumberSimpleCaseStudy class implements a simulation study around 
 * the exponential, square and log transformation of Gaussian variables.
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberSimpleCaseStudy {


	/**
	 * A population unit.<p>
	 * The variable x follows a uniform distribution (3,10).
	 * @author Mathieu Fortin - December 2023
	 */
	static class PopulationUnit {
		final double x;
		final double y;
		final int id;
		PopulationUnit(Matrix trueBeta, double trueStd, int id) {
			this.id = id;
			x = StatisticalUtility.getRandom().nextDouble() * 7 + 3;
			y = trueBeta.getValueAt(0, 0) + x * trueBeta.getValueAt(1, 0) + StatisticalUtility.getRandom().nextGaussian() * trueStd;
		}
	}

	/**
	 * A sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	@SuppressWarnings("serial")
	static class Sample extends ArrayList<PopulationUnit> {
		
		static Sample createSample(Matrix trueBeta, double trueVariance, int sampleSize) {
			Sample s = new Sample();
			double trueStd = Math.sqrt(trueVariance);
			for (int i = 0; i < sampleSize; i++) 
				s.add(new PopulationUnit(trueBeta, trueStd, i));
			return s;
		}
		
		Matrix getMatrixX() {
			Matrix xMat = new Matrix(size(), 2);
			for (int i = 0; i < size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, get(i).x);
			}
			return xMat;
		}

		Matrix getVectorY() {
			Matrix yVec = new Matrix(size(), 1);
			for (int i = 0; i < size(); i++) {
				yVec.setValueAt(i, 0, get(i).y);
			}
			return yVec;
		}
		
		Model getModel(Double sigma2) {
			Matrix xMat = getMatrixX();
			Matrix yVec = getVectorY();
			SymmetricMatrix invXtX = SymmetricMatrix.convertToSymmetricIfPossible(xMat.transpose().multiply(xMat).getInverseMatrix());
			Matrix betaHat = invXtX.multiply(xMat.transpose()).multiply(yVec);
			Matrix res = yVec.subtract(xMat.multiply(betaHat));
			int upsilon = yVec.m_iRows - xMat.m_iCols;
			double sigma2Hat = sigma2 != null ?
					sigma2 :
						res.transpose().multiply(res).getValueAt(0, 0) / upsilon;
			return new Model(this, betaHat, invXtX, sigma2Hat, upsilon);
		}
		
	}
	
	/**
	 * A model fitted to a sample of population units.
	 * @author Mathieu Fortin - December 2023
	 */
	static class Model {
		final Matrix betaHat;
		final SymmetricMatrix invXtX;
		final double sigma2Hat;
		final double sigmaHat;
		Matrix invXtXChol;
		final int upsilon;
		final Sample sample;
		final Matrix residuals;
		final Map<Integer, List<Integer>> sampleIndexMap;
		
		Model(Sample sample, Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.sample = sample;
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
			this.sigmaHat = Math.sqrt(sigma2Hat);
			residuals = sample.getVectorY().subtract(sample.getMatrixX().multiply(betaHat));
			sampleIndexMap = new HashMap<Integer, List<Integer>>();
		}
		
		private Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
		Matrix getRandomDeviate(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			Matrix betaHat_b = betaHat.add(betaDeviates);
			Matrix epsilon_b = StatisticalUtility.drawRandomVector(xMat.m_iRows, Type.GAUSSIAN).scalarMultiply(sigmaHat);
			Matrix result = xMat.multiply(betaHat_b).add(epsilon_b);
			return result.expMatrix();
		}
		
		Matrix createMatrixX(List<Double> xValues) {
			Matrix xMat = new Matrix(xValues.size(),2);
			for (int i = 0; i < xValues.size(); i++) {
				xMat.setValueAt(i, 0, 1d);
				xMat.setValueAt(i, 1, xValues.get(i));
			}
			return xMat;
		}
		
		Matrix getXBeta(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			return xMat.multiply(betaHat);
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

		private Matrix getBeauchampAndOlsonEstimatorVariance(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix xBeta = getXBeta(xValues);
			
			Matrix output = new Matrix(xValues.size(), 1);
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				double term1 = Math.exp(2*xBeta.getValueAt(ii, 0) + sigma2Hat);
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				int n = upsilon + xMat_ii.m_iCols;
				Matrix xInvXtXxT = xMat_ii.multiply(invXtX).multiply(xMat_ii.transpose());
				double phi = xInvXtXxT.getValueAt(0, 0) * n;
				double sigma2hat2 = sigma2Hat * sigma2Hat;
				double term2 = Math.exp(sigma2Hat) * (1 - 2 * sigma2Hat * (phi + 2*sigma2Hat)/n + 
						2* sigma2hat2 * (4*sigma2hat2 + 2 * (16d/3 + 2 * phi) * sigma2Hat + phi * phi + 4 * phi) / (n * n));
				double term3 = 1 - sigma2Hat/n * (sigma2Hat + 2 * phi) + sigma2hat2 / (2*n*n) * (sigma2hat2 + (16d/3 + 4*phi) * sigma2Hat + 4 * phi * phi + 8 * phi);
				output.setValueAt(ii, 0, term1 * (term2 - term3));
			}
			return output;
		}

		
		private Matrix getBaskervilleEstimator(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues);
			return xBeta.scalarAdd(sigma2Hat * .5).expMatrix();
		}
		
		private Matrix getBaskervilleEstimatorVariance(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues).scalarMultiply(2d).scalarAdd(sigma2Hat).expMatrix();
			return xBeta.scalarMultiply(Math.exp(sigma2Hat) - 1);
		}

		ComplexMatrix getComplexRandomDeviate(List<Double> xValues) {
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate - upsilon) / upsilon);
			ComplexNumber sigmaHat_b = sigma2Hat_b.sqrt();

			Matrix xMat = createMatrixX(xValues);

			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN));
			
			ComplexNumber[] cnArray = new ComplexNumber[xValues.size()];
			for (int ii = 0; ii < xValues.size(); ii++) {
				Matrix xMat_ii = xMat.getSubMatrix(ii, ii, 0, xMat.m_iCols - 1);
				ComplexNumber betaDeviates = new ComplexNumber(0, xMat_ii.multiply(betaDeviatesMat).getValueAt(0, 0));
				betaDeviates = betaDeviates.multiply(sigmaHat_b);
				Number epsilon_b = sigmaHat_b.multiply(StatisticalUtility.getRandom().nextGaussian());
				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
				ComplexNumber mean = betaDeviates.add(xBetaHat);
				ComplexNumber result = mean.add(epsilon_b);
				cnArray[ii] = result.exp();
			}
			return new ComplexMatrix(cnArray);
		}

	}
	
	private static void doRun(int sampleSize, int nbRealizations, double b0, double b1, double s2) throws IOException {
		System.out.println("Simulating [" + b0 + "; " + b1 + "; " + s2 +"] with sample size n = " + sampleSize);
		int nbInnerReal = 10000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberSimpleCaseStudy.class) + "caseStudy_" + s2 + "_" + sampleSize +".csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("x"));
		fields.add(new CSVField("MC_OLD"));
		fields.add(new CSVField("MC_OLDVar"));
		fields.add(new CSVField("MC_NEW"));
		fields.add(new CSVField("MC_NEW_imag"));
		fields.add(new CSVField("MC_NEWVar_real"));
		fields.add(new CSVField("MC_NEWVar_imag"));
		fields.add(new CSVField("Beauchamp"));
		fields.add(new CSVField("BeauchampVar"));
		fields.add(new CSVField("Baskerville"));
		fields.add(new CSVField("BaskervilleVar"));
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		trueBeta.setValueAt(0, 0, b0);
		trueBeta.setValueAt(1, 0, b1);
		double trueVariance = s2;
		
		List<Double> xValues = new ArrayList<Double>();
		for (double i = 3; i <= 10; i++) {
			xValues.add(i);
		}
		for (int real = 1; real <= nbRealizations; real++) {
			if (real % 1000 == 0) {
				System.out.println("Running realization " + real);
			}
			Sample s = Sample.createSample(trueBeta, trueVariance, sampleSize);
			Model m = s.getModel(null); // null means we are relying on the estimated variance
			Matrix beauchampAndOlsonEstimator = m.getBeauchampAndOlsonEstimator(xValues);
			Matrix beauchampAndOlsonEstimatorVar = m.getBeauchampAndOlsonEstimatorVariance(xValues);
			Matrix baskervilleEstimator = m.getBaskervilleEstimator(xValues);
			Matrix baskervilleEstimatorVar = m.getBaskervilleEstimatorVariance(xValues);
						
			MonteCarloEstimate mcEstimator = new MonteCarloEstimate();
			ComplexMonteCarloEstimate cmcEstimator = new ComplexMonteCarloEstimate();

			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				mcEstimator.addRealization(m.getRandomDeviate(xValues)); 
				cmcEstimator.addRealization(m.getComplexRandomDeviate(xValues)); 
			}
			
			Matrix mcMean = mcEstimator.getMean();
			Matrix mcVar = mcEstimator.getVariance();
			
			ComplexMatrix cmcMean = cmcEstimator.getMean();
//			Matrix cmcRealPart = cmcEstimator.getVarianceRealPart();
//			Matrix cmcImagPart = cmcEstimator.getVarianceImaginaryPart();

			ComplexSymmetricMatrix cmPseudoVariance = cmcEstimator.getPseudoVariance();
			
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				Object[] record = new Object[fields.size()];
				record[0] = real;
				record[1] = m.betaHat.getValueAt(0, 0);
				record[2] = m.betaHat.getValueAt(1, 0);
				record[3] = m.sigma2Hat;
				record[4] = xValues.get(ii);
				record[5] = mcMean.getValueAt(ii, 0);
				record[6] = mcVar.getValueAt(ii, ii);
				record[7] = cmcMean.getValueAt(ii, 0).realPart;
				record[8] = cmcMean.getValueAt(ii, 0).imaginaryPart;
//				record[9] = cmcRealPart.getValueAt(ii, 0);
//				record[10] = cmcImagPart.getValueAt(ii, 0);
				record[9] = cmPseudoVariance.getValueAt(ii, ii).realPart;
				record[10] = cmPseudoVariance.getValueAt(ii, ii).imaginaryPart;
				record[11] = beauchampAndOlsonEstimator.getValueAt(ii, 0);
				record[12] = beauchampAndOlsonEstimatorVar.getValueAt(ii, 0);
				record[13] = baskervilleEstimator.getValueAt(ii, 0);
				record[14] = baskervilleEstimatorVar.getValueAt(ii, 0);
				writer.addRecord(record);
			}
		}
		writer.close();
	}
	

	public static void main(String[] arg) throws IOException {
		doRun(25, 1000, 2, 0.25, 2);
//		doRun(Transformation.Exp, 50, 10000);
//		doRun(Transformation.Exp, 100, 10000);
//		doRun(Transformation.Exp, 200, 10000);
//		doRun(Transformation.Exp, 100);
//		doRun(Transformation.Exp, 200);
//		doRun(Transformation.Sqr, 50);
//		doRun(Transformation.Sqr, 100);
//		doRun(Transformation.Sqr, 200);
	}	
}
