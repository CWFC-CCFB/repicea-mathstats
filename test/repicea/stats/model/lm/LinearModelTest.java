/*
 * This file is part of the repicea-mathstats library.
 *
 * Copyright (C) 2023 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin - Canadian Wood Fibre Centre, CFS
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
package repicea.stats.model.lm;

import java.security.InvalidParameterException;

import org.junit.Assert;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.GaussianEstimate;
import repicea.stats.estimates.MonteCarloEstimate;
import repicea.stats.estimates.TruncatedGaussianEstimate;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.lm.LogBackTransformation.Estimator;
import repicea.util.ObjectUtility;

public class LinearModelTest {

	@Test
	public void simpleLinearModelWithOLSEstimator() throws Exception {
		String filename = ObjectUtility.getPackagePath(getClass()) + "datasetSingleObs.csv";
		DataSet ds = new DataSet(filename, true);
		LinearModel lm = new LinearModel(ds, "yTrans ~ lnDt_corr + BAL + dbhCm");
		lm.doEstimation();
		System.out.println(lm.getSummary());
		double expectedIntercept = -0.7229070890685245;
		Assert.assertEquals("Testing intercept estimate", expectedIntercept, lm.getParameters().getValueAt(0, 0), 1E-8);
		double expectedInterceptStd = 0.05541272279088098;
		double actualInterceptStd = Math.sqrt(lm.getEstimator().getParameterEstimates().getVariance().getValueAt(0, 0));
		Assert.assertEquals("Testing intercept estimate std", expectedInterceptStd, actualInterceptStd, 1E-8);
		double expectedPrediction = 0.8096288447607224;
		double actualPrediction = lm.getPredicted().getValueAt(0, 0);
		Assert.assertEquals("Testing first prediction", expectedPrediction, actualPrediction, 1E-8);
		Matrix newData = new Matrix(1,4);
		newData.setValueAt(0, 0, 1d);
		newData.setValueAt(0, 1, Math.log(5));
		newData.setValueAt(0, 2, 10d);
		newData.setValueAt(0, 3, 15d);
		double expectedNewPrediction = 0.33750837874909523;
		double actualNewPrediction = lm.getPredicted(newData).getValueAt(0, 0);
		Assert.assertEquals("Testing new prediction", expectedNewPrediction, actualNewPrediction, 1E-8);
	}
	
	@Test
	public void simpleLinearModelWithMaximumLikelihoodEstimator() throws Exception {
		String filename = ObjectUtility.getPackagePath(getClass()) + "datasetSingleObs.csv";
		DataSet ds = new DataSet(filename, true);
		Matrix startingValues = new Matrix(5, 1);
		startingValues.setValueAt(0, 0, -.5);
		startingValues.setValueAt(1, 0, .5);
		startingValues.setValueAt(2, 0, -0.007);
		startingValues.setValueAt(3, 0, -0.007);
		startingValues.setValueAt(4, 0, 0.3);
		LinearModel lm = new LinearModel(ds, "yTrans ~ lnDt_corr + BAL + dbhCm", startingValues);
		Assert.assertTrue("Testing the estimator is maximum likelihood", lm.getEstimator() instanceof MaximumLikelihoodEstimator);
		lm.doEstimation();
		System.out.println(lm.getSummary());
		
		double sigma2 = lm.getResidualVariance();
		Assert.assertEquals("Testing residual variance", 0.22291524, sigma2, 1E-8);
		
		double expectedIntercept = -0.72290709;
		Assert.assertEquals("Testing intercept estimate", expectedIntercept, lm.getParameters().getValueAt(0, 0), 1E-8);
		double expectedInterceptStd = 0.05539291;
		double actualInterceptStd = Math.sqrt(lm.getEstimator().getParameterEstimates().getVariance().getValueAt(0, 0));
		Assert.assertEquals("Testing intercept estimate std", expectedInterceptStd, actualInterceptStd, 1E-8);
		
		Matrix pred = lm.getPredicted();
		GaussianEstimate estimate = new GaussianEstimate(pred.getValueAt(0, 0), lm.getResidualVariance());
		MonteCarloEstimate mcEstimate = new MonteCarloEstimate();
		for (int i = 0; i < 1000000; i++) {
			mcEstimate.addRealization(estimate.getRandomDeviate().expMatrix().scalarAdd(-1d));
		}
//		Matrix actualPredOrig = lm.getPredOnLogBackTransformedScale(1d, true);
		Matrix actualPredOrig = LogBackTransformation.getMeanPredictedValuesOnOriginalScale(lm, 1d, Estimator.Naive);
		double expectedMean = mcEstimate.getMean().getValueAt(0, 0);
		Assert.assertEquals("Testing mean on original scale", expectedMean, actualPredOrig.getValueAt(0, 0), 1E-2);
		Matrix actualResidualVariance = LogBackTransformation.getResidualVariancesOnOriginalScale(lm, Estimator.Naive);
		double expectedVariance = mcEstimate.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing variance on original scale", expectedVariance, actualResidualVariance.getValueAt(0, 0), 2E-2);
	}

	@Test
	public void linearModelWithMaximumLikelihoodEstimatorButInconsistentStartingValues() throws Exception {
		String filename = ObjectUtility.getPackagePath(getClass()) + "datasetSingleObs.csv";
		DataSet ds = new DataSet(filename, true);
		Matrix startingValues = new Matrix(4, 1);
		startingValues.setValueAt(0, 0, -.5);
		startingValues.setValueAt(1, 0, .5);
		startingValues.setValueAt(2, 0, -0.007);
		startingValues.setValueAt(3, 0, -0.007);
		try {
			new LinearModel(ds, "yTrans ~ lnDt_corr + BAL + dbhCm", startingValues);
		} catch (InvalidParameterException e) {
			System.out.println(e.getMessage());
			return;
		}
		Assert.fail("Should have thrown an invalid parameter exception since the starting parameters are not consistent!");
	}

	@Test
	public void truncatedLinearModel() throws Exception {
		String filename = ObjectUtility.getPackagePath(getClass()) + "datasetSingleObs.csv";
		DataSet ds = new DataSet(filename, true);
		Matrix startingValues = new Matrix(5, 1);
		
		startingValues.setValueAt(0, 0, -0.72);
		startingValues.setValueAt(1, 0, 0.73);
		startingValues.setValueAt(2, 0, -0.005);
		startingValues.setValueAt(3, 0, -0.005);
		startingValues.setValueAt(4, 0, 0.22);

		LinearModelWithTruncatedGaussianErrorTerm lm = new LinearModelWithTruncatedGaussianErrorTerm(ds, "yTrans ~ lnDt_corr + BAL + dbhCm", startingValues, 0);
		lm.doEstimation();
		System.out.println(lm.getSummary());
		
		Matrix pred = lm.getPredicted();
		System.out.println("Predicted log scale = " + pred.getValueAt(0, 0));
		System.out.println("Variance log scale = " + lm.getResidualVariance());
		Assert.assertEquals("Testing residual variances", 0.32382955, lm.getResidualVariance(), 1E-8);
		
		
		Matrix parms = lm.getParameters();
		parms = parms.getSubMatrix(0, parms.m_iRows - 2, 0, 0);
		double xBeta = lm.getMatrixX().multiply(parms).getValueAt(0, 0);
		TruncatedGaussianEstimate estimate = new TruncatedGaussianEstimate(xBeta, lm.getResidualVariance());
		estimate.setLowerBoundValue(new Matrix(1,1));
		
		MonteCarloEstimate mcEstimate = new MonteCarloEstimate();
		MonteCarloEstimate mcEstimateOrig = new MonteCarloEstimate();
		for (int i = 0; i < 1000000; i++) {
			Matrix randomDeviate = estimate.getRandomDeviate();
			mcEstimate.addRealization(randomDeviate);
			mcEstimateOrig.addRealization(randomDeviate.expMatrix().scalarAdd(-1d));
		}
		
		double expectedPred = mcEstimate.getMean().getValueAt(0, 0);
		Assert.assertEquals("Testing pred on tranformed scale", expectedPred, pred.getValueAt(0, 0), 1E-2);
		
//		Matrix predOriginalScale = lm.getPredOnLogBackTransformedScale(1d, true); // offset = 1 true : with variance
		Matrix predOriginalScale = LogBackTransformation.getMeanPredictedValuesOnOriginalScale(lm, 1d, Estimator.Naive); 
		double expectedPredOrig = mcEstimateOrig.getMean().getValueAt(0, 0);
		Assert.assertEquals("Testing pred on original scale", expectedPredOrig, predOriginalScale.getValueAt(0, 0), 1E-2);
		Matrix resVarianceOriginalScale = LogBackTransformation.getResidualVariancesOnOriginalScale(lm, Estimator.Naive); 
		double expectedVarianceOrig = mcEstimateOrig.getVariance().getValueAt(0, 0);
		Assert.assertEquals("Testing variance on original scale", expectedVarianceOrig, resVarianceOriginalScale.getValueAt(0, 0), 2E-2);
		
		double expectedIntercept = -1.82587351;
		Assert.assertEquals("Testing intercept estimate", expectedIntercept, lm.getParameters().getValueAt(0, 0), 1E-8);
		double expectedInterceptStd = 0.09970163;
		double actualInterceptStd = Math.sqrt(lm.getEstimator().getParameterEstimates().getVariance().getValueAt(0, 0));
		Assert.assertEquals("Testing intercept estimate std", expectedInterceptStd, actualInterceptStd, 1E-8);
	}

	@Test
	public void truncatedLinearModelTest2() throws Exception {
		String filename = ObjectUtility.getPackagePath(getClass()) + "diameterIncrementData.csv";
		DataSet ds = new DataSet(filename, true);
		LinearModel lm = new LinearModel(ds, "yTrans ~ lnDt_corr + BAL*dbhCm.x + G_othersMinusBAL + meanFrostDay + meanCMI");
		lm.doEstimation();
		System.out.println(lm.getSummary());

		Matrix sigma2Res = new Matrix(1,1,lm.getResidualVariance(),0d);
		Matrix startingValues = lm.getParameters().matrixStack(sigma2Res, true);
		LinearModelWithTruncatedGaussianErrorTerm lmt = new LinearModelWithTruncatedGaussianErrorTerm(ds, 
				"yTrans ~ lnDt_corr + BAL*dbhCm.x + G_othersMinusBAL + meanFrostDay + meanCMI", 
				startingValues, 0);
		lmt.doEstimation();
		System.out.println(lmt.getSummary());

		double expectedIntercept = 0.4631388180861423;
		Assert.assertEquals("Testing intercept estimate", expectedIntercept, lm.getParameters().getValueAt(0, 0), 1E-8);
		double expectedInterceptStd = 0.04381822451441379 ;
		double actualInterceptStd = Math.sqrt(lm.getEstimator().getParameterEstimates().getVariance().getValueAt(0, 0));
		Assert.assertEquals("Testing intercept estimate std", expectedInterceptStd, actualInterceptStd, 1E-8);
	}

}
