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
import repicea.stats.estimators.MaximumLikelihoodEstimator;
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
		double expectedIntercept = -0.72290709;
		Assert.assertEquals("Testing intercept estimate", expectedIntercept, lm.getParameters().getValueAt(0, 0), 1E-8);
		double expectedInterceptStd = 0.05539291;
		double actualInterceptStd = Math.sqrt(lm.getEstimator().getParameterEstimates().getVariance().getValueAt(0, 0));
		Assert.assertEquals("Testing intercept estimate std", expectedInterceptStd, actualInterceptStd, 1E-8);
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
			LinearModel lm = new LinearModel(ds, "yTrans ~ lnDt_corr + BAL + dbhCm", startingValues);
		} catch (InvalidParameterException e) {
			System.out.println(e.getMessage());
			return;
		}
		Assert.fail("Should have thrown an invalid parameter exception since the starting parameters are not consistent!");
	}

}