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

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.serial.xml.XmlDeserializer;
import repicea.stats.Distribution.Type;
import repicea.stats.StatisticalUtility;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.util.ObjectUtility;

/**
 * The ComplexNumberDiffEquationCaseStudy class implements a simulation study on
 * differential growth equation. 
 * @author Mathieu Fortin - Oct 2023
 */
public class ComplexNumberGrowthModelCaseStudy {

	
	static class Tree implements Comparable<Tree> {
		final double dbhCm;
		final int id;
		double balM2ha;
		
		Tree(int id, double dbhCm) {
			this.dbhCm = dbhCm;
			this.id = id;
		}

		Tree(int id, double dbhCm, double balM2Ha) {
			this.id = id;
			this.dbhCm = dbhCm;
			this.balM2ha = balM2Ha;
		}

		double getBasalAreaM2Ha() {
			return dbhCm * dbhCm * Math.PI * 0.000625; // /40000 * 10000 / 400 = 0.000625
		}

		@Override
		public int compareTo(Tree o) {
			if (o.dbhCm > dbhCm)
				return 1;
			if (o.dbhCm == dbhCm) {
				return 0;
			} else {
				return -1;
			}
		}
	}

	static class Stand {
		final List<Tree> trees;
		final int age;
		Stand(int age, int nbTrees) {
			this.age = age;
			this.trees = new ArrayList<Tree>();
			for (int i = 0; i < nbTrees; i++) {
				trees.add(new Tree(i, StatisticalUtility.getRandom().nextDouble() * age + 2.5));
			}
		}
		
		void updateBAL() {
			Collections.sort(trees);
			double balM2Ha = 0;
			for (int i = 0; i < trees.size(); i++) {
				trees.get(i).balM2ha = balM2Ha;
				balM2Ha += trees.get(i).getBasalAreaM2Ha();
			}
		}
		
		double getBasalAreaM2Ha() {
			double sumBAHa = 0;
			for (Tree t : trees) {
				sumBAHa += t.getBasalAreaM2Ha();
			}
			return sumBAHa;
		}
		
		Stand updateStand(Matrix incrementCm) {
			Stand s = new Stand(age + 5, 0);
			for (int i = 0; i < incrementCm.m_iRows; i++) {
				Tree t = trees.get(i);
				double increment = incrementCm.getValueAt(i, 0);
				s.trees.add(new Tree(t.id, t.dbhCm + increment));
			}
			return s;
		}

		
	}
	
	static class GrowthModel {
		final Matrix beta;
		final double sigma2;
		final double sigma;
		SymmetricMatrix invXtX;
		Matrix chol;
		int upsilon;
		
		GrowthModel(Matrix beta, double sigma2) {
			this.beta = beta.getDeepClone();
			this.sigma2 = sigma2;
			this.sigma = Math.sqrt(sigma2);
		}

		GrowthModel(Matrix beta, double sigma2, int upsilon, SymmetricMatrix invXtX) {
			this(beta, sigma2);
			this.upsilon = upsilon;
			this.invXtX = invXtX;
		}
		
		
		Matrix getDbhIncrement(List<Tree> trees, Matrix newBeta) {
			Matrix x = getXMatrix(trees);
			Matrix pred = x.multiply(newBeta).add(StatisticalUtility.drawRandomVector(x.m_iRows, Type.GAUSSIAN).scalarMultiply(sigma));
			return pred.expMatrix().scalarAdd(-1d);
		}
		
		private Matrix getChol() {
			if (chol == null) {
				chol = invXtX.getLowerCholTriangle();
			}
			return chol;
		}
		
		MonteCarloEstimate runSimulation(int nbRealizations, int nbSteps, Stand initStand) {
			boolean isTrueModel = invXtX == null;
			MonteCarloEstimate estimate = new MonteCarloEstimate();
			Matrix beta_b = null;
			for (int real = 0; real < nbRealizations; real++) {
				Matrix realization = new Matrix(nbSteps, 1);
				if (!isTrueModel) {
					beta_b = beta.add(getChol().multiply(StatisticalUtility.drawRandomVector(beta.m_iRows, Type.GAUSSIAN)));
				}
				Stand s = initStand;
				s.updateBAL();
				realization.setValueAt(0, 0, s.getBasalAreaM2Ha());
				for (int j = 1; j < nbSteps; j++) {
					s = s.updateStand(getDbhIncrement(s.trees, beta_b == null ? beta : beta_b));
					s.updateBAL();
					realization.setValueAt(j, 0, s.getBasalAreaM2Ha());
				}
				estimate.addRealization(realization);
			}
			return estimate;
		}

		
	}
	
	static List<Tree> produceSample(int n) {
		List<Tree> trees = new ArrayList<Tree>();
		for (int i = 0; i < n; i++) {
			double dbhCm = 5 + StatisticalUtility.getRandom().nextDouble() * 40;
			double balM2Ha = 5 + StatisticalUtility.getRandom().nextDouble() * 20;
			trees.add(new Tree(i, dbhCm, balM2Ha));
		}
		return trees;
	}
	
	static Matrix getXMatrix(List<Tree> trees) {
		Matrix x = new Matrix(trees.size(), 3, 1, 0);
		for (int i = 0; i < trees.size(); i++) {
			x.setValueAt(i, 1, trees.get(i).balM2ha);
			x.setValueAt(i, 2, trees.get(i).dbhCm);
		}
		return x;
	}
	
	public static GrowthModel fitGrowthModel(GrowthModel trueModel, int sampleSize) {
		List<Tree> trees = produceSample(sampleSize);
		Matrix dbhIncrementCm = trueModel.getDbhIncrement(trees, trueModel.beta);
		Matrix y = dbhIncrementCm.scalarAdd(1d).logMatrix();
		Matrix x = getXMatrix(trees);

		SymmetricMatrix invXtX = SymmetricMatrix.convertToSymmetricIfPossible(x.transpose().multiply(x).getInverseMatrix());
		Matrix betaHat = invXtX.multiply(x.transpose()).multiply(y);
		Matrix res = y.subtract(x.multiply(betaHat));
		int upsilon = x.m_iRows - x.m_iCols;
		double sigma2Hat = res.transpose().multiply(res).getValueAt(0, 0) / (upsilon);
		GrowthModel m = new GrowthModel(betaHat, sigma2Hat, upsilon, invXtX);
		return m;
	}
	
	
	public static void main(String[] arg) throws IOException {
		Matrix trueBeta = new Matrix(3,1);
		trueBeta.setValueAt(0, 0, 0.8);
		trueBeta.setValueAt(1, 0, -0.02);
		trueBeta.setValueAt(2, 0, -0.02);
		double trueSigma2 = 0.25;
		GrowthModel trueModel = new GrowthModel(trueBeta, trueSigma2);
				
		int maxSteps = 15;
		int nbRealizations = 1000;
		int nbInnerRealizations = 1000;
		String filename = ObjectUtility.getPackagePath(ComplexNumberGrowthModelCaseStudy.class).replace("bin/", "") + "initStand.zml";
//		Stand initStand = new Stand(5, 50);
//		XmlSerializer serializer = new XmlSerializer(filename);
//		serializer.writeObject(initStand);
		XmlDeserializer deserializer = new XmlDeserializer(filename);
		Stand initStand = (Stand) deserializer.readObject();
		MonteCarloEstimate truth = trueModel.runSimulation(10000, maxSteps, initStand);
		MonteCarloEstimate estimatorRealizations = new MonteCarloEstimate();
		for (int i = 1; i <= nbRealizations; i++) {
			if (i%100 == 0) {
				System.out.println("Realization " + i + "/" + nbRealizations);
			}
			GrowthModel m = fitGrowthModel(trueModel, 100);
			MonteCarloEstimate estimate = m.runSimulation(nbInnerRealizations, maxSteps, initStand);
			estimatorRealizations.addRealization(estimate.getMean());
		}
		System.out.println("Reference = " + truth.getMean().toString());
		System.out.println("Expectation = " + estimatorRealizations.getMean().toString());
		System.out.println("Variance estimator = " + estimatorRealizations.getVariance().diagonalVector().toString());
	}	
}
