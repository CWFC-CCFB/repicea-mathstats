/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm.measerr;

import java.util.logging.ConsoleHandler;
import java.util.logging.Level;

import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.math.optimizer.NewtonRaphsonOptimizer;
import repicea.stats.data.DataSet;
import repicea.stats.estimates.Estimate;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class SIMEXTest {

	
	@BeforeClass
	public static void doThisBefore() {
		Level l = Level.OFF;
		NewtonRaphsonOptimizer.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		ConsoleHandler ch = new ConsoleHandler();
		ch.setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).addHandler(ch);		
		SIMEXModel.OverrideVarianceForTest = true;
	}
	
	@Test
	public void simpleTest() throws Exception {
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample0.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific");
		glm.doEstimation();
		SIMEXModel s = new SIMEXModel(glm, "distanceToConspecific", "variance");
		s.setNumberOfBootstrapRealizations(200);
		s.doEstimation();
		Assert.assertTrue("Checking if successfully extrapolated", s.getEstimator().isConvergenceAchieved());
		Estimate<?> estimate = s.getEstimator().getParameterEstimates();
		Assert.assertEquals("Checking parm estimate 1", -0.1687, estimate.getMean().getValueAt(0, 0), 4E-3);
		Assert.assertEquals("Checking standard error parm estimate 1", 0.067178, Math.sqrt(estimate.getVariance().getValueAt(0, 0)), 6E-4);
		Assert.assertEquals("Checking parm 2", -0.05505, estimate.getMean().getValueAt(1, 0), 4E-4);
		Assert.assertEquals("Checking standard error parm estimate 2", 0.004364, Math.sqrt(estimate.getVariance().getValueAt(1, 1)), 6E-5);
	}
	
	@Test
	public void multipleOccurrenceOfMeasurementErrorTest() throws Exception {
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample0.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific + log(10 + distanceToConspecific)");
		glm.doEstimation();
		SIMEXModel s = new SIMEXModel(glm, "distanceToConspecific", "variance");
		s.setNumberOfBootstrapRealizations(200);
		s.doEstimation();
		Assert.assertTrue("Checking if successfully extrapolated", s.getEstimator().isConvergenceAchieved());
		Estimate<?> estimate = s.getEstimator().getParameterEstimates();
		Assert.assertEquals("Checking parm 2", -0.062640, estimate.getMean().getValueAt(1, 0), 5E-3);
		Assert.assertEquals("Checking standard error parm estimate 2", 0.020396, Math.sqrt(estimate.getVariance().getValueAt(1, 1)), 1E-3);
	}

	@AfterClass
	public static void doThatAfter() {
		SIMEXModel.OverrideVarianceForTest = false;
	}

	@Test
	public void checkingPredParmsDataSetTest() throws Exception {
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample0.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific");
		glm.doEstimation();
		SIMEXModel s = new SIMEXModel(glm, "distanceToConspecific", "variance");
		s.setNumberOfBootstrapRealizations(100);
		s.doEstimation();
		DataSet dsPred = s.getPredictedParameterEstimates();
		Assert.assertTrue("Testing that the dataset is not null", dsPred != null);
		int index = dsPred.getIndexOfThisField("pred");
		double parmValue = Double.parseDouble(dsPred.getObservations().get(0).getValueAt(index).toString());
		Estimate<?> estimate = s.getEstimator().getParameterEstimates();
		Assert.assertEquals("Checking parm estimate 1", estimate.getMean().getValueAt(0, 0), parmValue, 1E-8);
	}

	@Test
	public void checkingSIMEXPredictionsTest() throws Exception {
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample0.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific");
		glm.doEstimation();
		SIMEXModel s = new SIMEXModel(glm, "distanceToConspecific", "variance");
		s.setNumberOfBootstrapRealizations(100);
		s.doEstimation();
		Matrix pred = s.getPredicted();
		Assert.assertTrue("Testing that the dataset is not null", pred != null);
		double predValue = pred.getValueAt(0, 0);
		Assert.assertEquals("Checking parm estimate 1", 0.054082156776618984, predValue, 3E-3);
	}
	
	@Test
	public void checkingSIMEXWithNegativeBinomial() throws Exception {
		boolean before = SIMEXModel.OverrideVarianceForTest;
		SIMEXModel.OverrideVarianceForTest = false;
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("recruitEPR.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, 
				GLMDistribution.NegativeBinomial, 
				Type.Log, "y ~ TotalPrcp + G_F + G_R + occIndex10km + timeSince1970");
		glm.doEstimation();
		System.out.println(glm.getSummary());
		SIMEXModel s = new SIMEXModel(glm, "occIndex10km", "occIndex10kmVar");
		s.setNumberOfBootstrapRealizations(200);
		s.doEstimation();
		System.out.println(s.getSummary());
		Matrix pred = s.getPredicted();
		Assert.assertTrue("Testing that the dataset is not null", pred != null);
		double parmEst = s.getParameters().getValueAt(0, 0);
		Assert.assertEquals("Checking parm estimate: intercept", -2.86, parmEst, 2E-1);
		parmEst = s.getParameters().getValueAt(1, 0);
		Assert.assertEquals("Checking parm estimate: TotalPrcp", 0.00155, parmEst, 3E-4);
		parmEst = s.getParameters().getValueAt(2, 0);
		Assert.assertEquals("Checking parm estimate: G_F", -0.1165, parmEst, 3E-3);
		parmEst = s.getParameters().getValueAt(3, 0);
		Assert.assertEquals("Checking parm estimate: G_R", -0.0713, parmEst, 2E-3);
		parmEst = s.getParameters().getValueAt(4, 0);
		Assert.assertEquals("Checking parm estimate: occIndex10km", 3.37, parmEst, 2E-1);
		parmEst = s.getParameters().getValueAt(5, 0);
		Assert.assertEquals("Checking parm estimate: timeSince1970", 0.0450, parmEst, 2E-3);
		parmEst = s.getParameters().getValueAt(6, 0);
		Assert.assertEquals("Checking parm estimate: theta", 1.83, parmEst, 1E-1);
		SIMEXModel.OverrideVarianceForTest = before; 
	}


	
}
