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
import java.util.List;

import repicea.io.FormatField;
import repicea.io.javacsv.CSVField;
import repicea.io.javacsv.CSVWriter;
import repicea.stats.Distribution;
import repicea.stats.StatisticalUtility;
import repicea.util.ObjectUtility;

public class ComplexNumberCaseStudy2 {

	
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
			SymmetricMatrix omega = invXtX.scalarMultiply(sigma2Hat);
			return new Model(betaHat, invXtX, sigma2Hat, upsilon);
		}
		
	}
	
	static class Model {
		final Matrix betaHat;
		final SymmetricMatrix invXtX;
		final double sigma2Hat;
		final double sigmaHat;
		Matrix invXtXChol;
		final int upsilon;
		
		Model(Matrix betaHat, SymmetricMatrix invXtX, double sigma2Hat, int upsilon) {
			this.betaHat = betaHat;
			this.invXtX = invXtX; 
			this.sigma2Hat = sigma2Hat;
			this.upsilon = upsilon;
			this.sigmaHat = Math.sqrt(sigma2Hat);
		}
		
		private Matrix getOmegaChol() {
			if (invXtXChol == null) {
				invXtXChol = invXtX.getLowerCholTriangle();
			}
			return invXtXChol;
		}
		
		double getRandomDeviate(double x) {
			Matrix xMat = new Matrix(1,2);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, x);
			Matrix betaDeviates = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN)).scalarMultiply(sigmaHat);
			
			Matrix betaHat_b = betaHat.add(betaDeviates);
			double epsilon_b = StatisticalUtility.getRandom().nextGaussian() * sigmaHat;
			return betaHat_b.getValueAt(0, 0) + betaHat_b.getValueAt(1, 0) * x + epsilon_b;
		}
		
		ComplexNumber getComplexRandomDeviate(double x) {
			Matrix xMat = new Matrix(1,2);
			xMat.setValueAt(0, 0, 1d);
			xMat.setValueAt(0, 1, x);
			
			double chiSquareDeviate = StatisticalUtility.getRandom().nextChiSquare(upsilon);
			ComplexNumber sigma2Hat_b = new ComplexNumber(sigma2Hat, sigma2Hat * (chiSquareDeviate - upsilon) / upsilon);
			Number sigmaHat_b = sigma2Hat_b.sqrt();
			
			Matrix betaDeviatesMat = getOmegaChol().multiply(StatisticalUtility.drawRandomVector(betaHat.m_iRows, Distribution.Type.GAUSSIAN));
			ComplexNumber betaDeviates = new ComplexNumber(0, xMat.multiply(betaDeviatesMat).getValueAt(0, 0));
			betaDeviates = betaDeviates.multiply(sigmaHat_b);
			
			Number epsilon_b = sigmaHat_b instanceof ComplexNumber ? 
					((ComplexNumber) sigmaHat_b).multiply(StatisticalUtility.getRandom().nextGaussian()) :
						sigmaHat_b.doubleValue() * StatisticalUtility.getRandom().nextGaussian();
			
			Number mean = betaDeviates.add(betaHat.getValueAt(0, 0) + betaHat.getValueAt(1, 0) * x);
			return (ComplexNumber) ((ComplexNumber) mean).add(epsilon_b);
		}

	}
	

	public static void main(String[] arg) throws IOException {
		int nbRealizations = 50000;
		int nbInnerReal = 10000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberCaseStudy2.class) + "caseStudy2.csv";
		filename = filename.replace("bin/", "");
		CSVWriter writer = new CSVWriter(new File(filename), false);
		List<FormatField> fields = new ArrayList<FormatField>();
		fields.add(new CSVField("RealID"));
		fields.add(new CSVField("b0"));
		fields.add(new CSVField("b1"));
		fields.add(new CSVField("sigma2"));
		fields.add(new CSVField("y2OLD"));
		fields.add(new CSVField("y2NEW"));
		writer.setFields(fields);
		
		Matrix trueBeta = new Matrix(2,1);
		trueBeta.setValueAt(0, 0, 2);
		trueBeta.setValueAt(1, 0, 1.5);
		double trueVariance = 4d;
		
		for (int real = 0; real < nbRealizations; real++) {
			System.out.println("Running realization " + (real + 1));
			Sample s = Sample.createSample(trueBeta, trueVariance, 30);
			Model m = s.getModel(null); // null means we are relying on the estimated variance
			double sumOLD = 0;
			double sumNEW = 0;
			for (int innerReal = 0; innerReal < nbInnerReal; innerReal++) {
				Number n1 = m.getRandomDeviate(5);
				sumOLD += Math.exp(n1.doubleValue());
				ComplexNumber n2 = m.getComplexRandomDeviate(5);
				sumNEW += n2.exp().doubleValue();
			}
			sumOLD /= nbInnerReal;
			sumNEW /= nbInnerReal;
			Object[] record = new Object[6];
			record[0] = real;
			record[1] = m.betaHat.getValueAt(0, 0);
			record[2] = m.betaHat.getValueAt(1, 0);
			record[3] = m.sigma2Hat;
			record[4] = sumOLD;
			record[5] = sumNEW;
			writer.addRecord(record);
		}
		writer.close();
	}	
}
