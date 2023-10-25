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
import java.security.InvalidParameterException;
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
import repicea.stats.sampling.SamplingUtility;
import repicea.util.ObjectUtility;

/**
 * The ComplexNumberSimpleCaseStudy class implements a simulation study around 
 * the exponential, square and log transformation of Gaussian variables.
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberSimpleCaseStudy {

	
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
		
		private Matrix sampleResiduals(int n) {
			if (!sampleIndexMap.containsKey(n)) {
				List<Integer> sampleIndex = new ArrayList<Integer>();
				for (int i = 0; i < n; i++)
					sampleIndex.add(i);
				sampleIndexMap.put(n, sampleIndex);
			}
			
			List<Integer> index = SamplingUtility.getSample(sampleIndexMap.get(n), n, true);
			Matrix sampledResiduals = residuals.getSubMatrix(index, null, false); // we do not sort the indices
			return sampledResiduals;
		}
		
		Matrix getRandomDeviate(List<Double> xValues, Transformation t, boolean useNonParametricMethod) {
			Matrix xMat = createMatrixX(xValues);
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			Matrix betaHat_b = betaHat.add(betaDeviates);
			
			Matrix epsilon_b = useNonParametricMethod ? 
					sampleResiduals(xMat.m_iRows) :
						StatisticalUtility.drawRandomVector(xMat.m_iRows, Type.GAUSSIAN).scalarMultiply(sigmaHat);
			
			Matrix result = xMat.multiply(betaHat_b).add(epsilon_b);
			return convertPredictionToOriginalScale(result, t);
		}
		
		private Matrix convertPredictionToOriginalScale(Matrix transScalePred, Transformation t) {
			switch(t) {
			case Sqr:
				return transScalePred.elementWisePower(2d);
			case Exp:
				return transScalePred.expMatrix();
			case Log:
				return transScalePred.logMatrix();
			case Sqrt:
				return transScalePred.elementWisePower(.5);
			default:
				throw new InvalidParameterException("This transformation " + t.name() + " is not supported!");
			}
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

		
		private Matrix getBaskervilleEstimator(List<Double> xValues) {
			Matrix xBeta = getXBeta(xValues);
			return xBeta.scalarAdd(sigma2Hat * .5).expMatrix();
		}
		
		private Matrix getSmearingEstimate(List<Double> xValues, Transformation t) {
			Matrix xBeta = getXBeta(xValues);
			Matrix res_hat = sample.getVectorY().subtract(sample.getMatrixX().multiply(betaHat));
			Matrix result = null;
			for (int i = 0; i < res_hat.m_iRows; i++) {
				result = result == null ?
						convertPredictionToOriginalScale(xBeta.scalarAdd(res_hat.getValueAt(i, 0)), t) :
							result.add(convertPredictionToOriginalScale(xBeta.scalarAdd(res_hat.getValueAt(i, 0)), t));	
			}
			return result.scalarMultiply(1d / res_hat.m_iRows);
		}
		
		ComplexNumber[] getComplexRandomDeviate(List<Double> xValues, Transformation t, boolean useNonParametricMethod) {
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

				Number epsilon_b = useNonParametricMethod ?
						sampleResiduals(1).getValueAt(0, 0) :
							sigmaHat_b.multiply(StatisticalUtility.getRandom().nextGaussian());
				double xBetaHat = xMat_ii.multiply(betaHat).getValueAt(0, 0); 
				ComplexNumber mean = betaDeviates.add(xBetaHat);
				ComplexNumber result = mean.add(epsilon_b);
				switch(t) {
				case Sqr:
					cnArray[ii] = result.square();
					break;
				case Exp:
					cnArray[ii] = result.exp();
					break;
				case Log:
					cnArray[ii] = result.log();
					break;
				case Sqrt:
					cnArray[ii] = result.sqrt();
					break;
				default:
					throw new InvalidParameterException("The transformation " + t.name() + " is not supported!");
				}
			}
			return cnArray;
		}

		private Matrix getGregoireEstimator(List<Double> xValues) {
			Matrix xMat = createMatrixX(xValues);
			Matrix varXBeta = xMat.multiply(invXtX).multiply(xMat.transpose()).scalarMultiply(sigma2Hat).diagonalVector();
			Matrix pred = getXBeta(xValues);
			
			return pred.elementWisePower(2d).subtract(varXBeta).scalarAdd(sigma2Hat);
		}
	}
	
	static enum Transformation {Sqr, Exp, Log, Sqrt}
	
	private static void addComplexNumberArray(ComplexNumber[] originalArray, ComplexNumber[] newArray) {
		for (int ii = 0; ii < originalArray.length; ii++) {
			originalArray[ii] = originalArray[ii].add(newArray[ii]);
		}
	}
	

	private static void doRun(Transformation t, int sampleSize) throws IOException {
		System.out.println("Simulating " + t.name() + " with sample size n = " + sampleSize);
		int nbRealizations = 10000;
		int nbInnerReal = 10000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberSimpleCaseStudy.class) + "caseStudy_" + t.name() + sampleSize +".csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("x"));
		fields.add(new CSVField("MC_OLD"));
//		fields.add(new CSVField("MC_OLD_np"));
		fields.add(new CSVField("MC_NEW"));
		fields.add(new CSVField("MC_NEW_imag"));
//		fields.add(new CSVField("MC_NEW_np"));
//		fields.add(new CSVField("MC_NEW_np_imag"));
		if (t == Transformation.Exp) {
			fields.add(new CSVField("Beauchamp"));
			fields.add(new CSVField("Baskerville"));
			fields.add(new CSVField("Smearing"));
		} else if (t == Transformation.Sqr) {
			fields.add(new CSVField("Smearing"));
			fields.add(new CSVField("Gregoire"));
		}
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		double trueVariance;
		if (t == Transformation.Exp) {
			trueBeta.setValueAt(0, 0, 2);
			trueBeta.setValueAt(1, 0, 0.25);
			trueVariance = 1d;
		} else if (t == Transformation.Sqr) {  
			trueBeta.setValueAt(0, 0, 2);
			trueBeta.setValueAt(1, 0, 0.25);
			trueVariance = 2d;
		} else {
			trueBeta.setValueAt(0, 0, 2.5);
			trueBeta.setValueAt(1, 0, 4);
			trueVariance = 12d;
		}
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
			Matrix beauchampAndOlsonEstimator = null;
			Matrix baskervilleEstimator = null;
			Matrix smearingEstimate = null;
			Matrix gregoireEstimator = null;
			if (t == Transformation.Exp) {
				beauchampAndOlsonEstimator = m.getBeauchampAndOlsonEstimator(xValues);
				baskervilleEstimator = m.getBaskervilleEstimator(xValues);
				smearingEstimate = m.getSmearingEstimate(xValues, t);
			} else if (t == Transformation.Sqr) {
				smearingEstimate = m.getSmearingEstimate(xValues, t);
				gregoireEstimator = m.getGregoireEstimator(xValues);
			}
						
			Matrix sumOLD = new Matrix(xValues.size(), 1);
//			Matrix sumOLD_np = new Matrix(xValues.size(), 1);
			ComplexNumber[] sumNEW = new ComplexNumber[xValues.size()];
			for (int ii = 0; ii < sumNEW.length; ii++) {
				sumNEW[ii] = new ComplexNumber(0,0);
			}
//			ComplexNumber[] sumNEW_np = new ComplexNumber[xValues.size()];
//			for (int ii = 0; ii < sumNEW.length; ii++) {
//				sumNEW_np[ii] = new ComplexNumber(0,0);
//			}

			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				sumOLD = sumOLD.add(m.getRandomDeviate(xValues, t, false)); // false: regular parametric method
//				sumOLD_np = sumOLD_np.add(m.getRandomDeviate(xValues, t, true)); // true: non parametric method
				ComplexNumber[] n2 = m.getComplexRandomDeviate(xValues, t, false); // false : regular parametric method
				addComplexNumberArray(sumNEW, n2);
//				ComplexNumber[] n2_np = m.getComplexRandomDeviate(xValues, t, true); // true : non parametric method
//				addComplexNumberArray(sumNEW_np, n2_np);
			}
			
			sumOLD = sumOLD.scalarMultiply(1d / nbInnerReal);
//			sumOLD_np = sumOLD_np.scalarMultiply(1d / nbInnerReal);
			for (int ii = 0; ii < sumNEW.length; ii++) {
				sumNEW[ii] = sumNEW[ii].multiply(1d / nbInnerReal);
			}
//			for (int ii = 0; ii < sumNEW_np.length; ii++) {
//				sumNEW_np[ii] = sumNEW_np[ii].multiply(1d / nbInnerReal);
//			}
			
			for (int ii = 0; ii < xValues.size(); ii++) {
				Object[] record = new Object[fields.size()];
				record[0] = real;
				record[1] = m.betaHat.getValueAt(0, 0);
				record[2] = m.betaHat.getValueAt(1, 0);
				record[3] = m.sigma2Hat;
				record[4] = xValues.get(ii);
				record[5] = sumOLD.getValueAt(ii, 0);
//				record[6] = sumOLD_np.getValueAt(ii, 0);
				record[6] = sumNEW[ii].doubleValue();
				record[7] = sumNEW[ii].imaginaryPart;
//				record[9] = sumNEW_np[ii].doubleValue();
//				record[10] = sumNEW_np[ii].imaginaryPart;
				if (t == Transformation.Exp) {
					record[8] = beauchampAndOlsonEstimator.getValueAt(ii, 0);
					record[9] = baskervilleEstimator.getValueAt(ii, 0);
					record[10] = smearingEstimate.getValueAt(ii, 0);
				} else  if (t == Transformation.Sqr) {
					record[8] = smearingEstimate.getValueAt(ii, 0);
					record[9] = gregoireEstimator.getValueAt(ii, 0);
				}
				writer.addRecord(record);
			}
		}
		writer.close();

	}
	
	public static void main(String[] arg) throws IOException {
		doRun(Transformation.Exp, 50);
		doRun(Transformation.Exp, 100);
		doRun(Transformation.Exp, 200);
		doRun(Transformation.Sqr, 50);
		doRun(Transformation.Sqr, 100);
		doRun(Transformation.Sqr, 200);
//		doRun(Transformation.Log, 50);
//		doRun(Transformation.Log, 100);
//		doRun(Transformation.Log, 200);
//		doRun(Transformation.Sqrt, 50);
//		doRun(Transformation.Sqrt, 100);
//		doRun(Transformation.Sqrt, 200);
	}	
}
