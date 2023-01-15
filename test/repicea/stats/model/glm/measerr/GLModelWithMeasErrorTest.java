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

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import repicea.math.optimizer.AbstractOptimizer.LineSearchMethod;
import repicea.math.optimizer.NewtonRaphsonOptimizer;
import repicea.stats.data.DataSet;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.GeneralizedLinearModel;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class GLModelWithMeasErrorTest {

	@BeforeClass
	public static void doThis() {
		Level l = Level.OFF;
		NewtonRaphsonOptimizer.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		ConsoleHandler ch = new ConsoleHandler();
		ch.setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).addHandler(ch);		
	}

//	@Test
//	public void testAlteredGaussianFunctionValue() {
//		double sigma2 = 3;
//		AlteredGaussianFunction afg = new AlteredGaussianFunction(sigma2);
//		double mu = 2;
//		afg.setVariableValue(afg.X_INDEX, mu);
//		double w = 5;
//		afg.wValue = w;
//		double observedValue = afg.getValue();
//		double expectedValue = GaussianUtility.getProbabilityDensity(w, mu, sigma2);
//		Assert.assertEquals("Comparing densities", expectedValue, observedValue, 1E-8);
//	}
//
//	@Test
//	public void testAlteredGaussianFunctionGradient() {
//		double sigma2 = 3;
//		AlteredGaussianFunction afg = new AlteredGaussianFunction(sigma2);
//		double mu = 2;
//		afg.setVariableValue(afg.X_INDEX, mu);
//		double w = 5;
//		afg.wValue = w;
//		double observedValue = afg.getGradient().getValueAt(0, 0);
//		double df_dsigma2 = ((w - mu) * (w - mu) - sigma2) * Math.exp(-(w - mu) * (w - mu) / (2 * sigma2)) / (2 * Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 2.5));
//
//		Assert.assertEquals("Comparing gradients", df_dsigma2, observedValue, 1E-8);
//	}
//	
//	@Test
//	public void testAlteredGaussianFunctionHessian() {
//		double sigma2 = 3;
//		AlteredGaussianFunction afg = new AlteredGaussianFunction(sigma2);
//		double mu = 2;
//		afg.setVariableValue(afg.X_INDEX, mu);
//		double w = 5;
//		afg.wValue = w;
//		double observedValue = afg.getHessian().getValueAt(0, 0);
//		double d2f_d2sigma2 = (3 * sigma2 * sigma2 - 6 * sigma2 * (w - mu) * (w -mu) + (w - mu) * (w -mu) * (w - mu) * (w -mu)) *
//				Math.exp(-(w - mu) * (w - mu) / (2 * sigma2)) / (4 * Math.sqrt(2 * Math.PI) * Math.pow(sigma2, 4.5));
//
//		Assert.assertEquals("Comparing hessians", d2f_d2sigma2, observedValue, 1E-8);
//	}

	
	
//	@Ignore
//	@Test
//    public void TestGLModelWithMeasurementError() throws Exception {
// 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("OccurrencePartDataset_ERS.csv");
//		DataSet dataSet = new DataSet(filename, true);
//		
//		GLMWithMeasurementError glmWithMeasError = new GLMWithMeasurementError(dataSet, "occurred ~ lnDt + distanceToConspecificOLD",
//				new GLMUniformBerksonMeasErrorDefinition("distanceToConspecificOLD", 0d, null, null, .1));
//		((MaximumLikelihoodEstimator) glmWithMeasError.getEstimator()).setLineSearchMethod(LineSearchMethod.HALF_STEP);
//		glmWithMeasError.doEstimation();
//		glmWithMeasError.getSummary();
//		Assert.assertTrue("Testing if convergence has been reached", glmWithMeasError.getEstimator().isConvergenceAchieved());
//		Assert.assertEquals("Testing the parameter estimate", -13.911668422111854, glmWithMeasError.getParameters().getValueAt(2, 0), 1E-8);
//	}

	@Test
    public void TestGLModelWithClassicalMeasurementError() throws Exception {
 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample0.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.CLogLog, "y ~ distanceToConspecific");
		glm.doEstimation();
		glm.getSummary();
		GLMWithMeasurementError glmWithMeasError = new GLMWithMeasurementError(dataSet, "y ~ distanceToConspecific",
				glm.getParameters(), new GLMNormalClassicalMeasErrorDefinition("distanceToConspecific", "variance", .2));
		((MaximumLikelihoodEstimator) glmWithMeasError.getEstimator()).setLineSearchMethod(LineSearchMethod.HALF_STEP);
		glmWithMeasError.doEstimation();
		glmWithMeasError.getSummary();
		Assert.assertTrue("Testing if convergence has been reached", glmWithMeasError.getEstimator().isConvergenceAchieved());
		Assert.assertEquals("Testing the parameter estimate", -0.054727557137678795, glmWithMeasError.getParameters().getValueAt(1, 0), 1E-8);
	}


//	@Test
//    public void TestGLModelWithClassicalMeasurementError2() throws Exception {
// 		String filename = ObjectUtility.getPackagePath(GLModelWithMeasErrorTest.class).concat("sample1.csv");
//		DataSet dataSet = new DataSet(filename, true);
//		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, Type.CLogLog, "y ~ trueDistanceToConspecific");
//		glm.doEstimation();
//		glm.getSummary();
//		Matrix initParms = new Matrix(2,1);
//		initParms.setValueAt(0, 0, .2);
//		initParms.setValueAt(1, 0, -.2);
//		GLMWithMeasurementError glmWithMeasError = new GLMWithMeasurementError(dataSet, "y ~ distanceToConspecific",
//				initParms, new GLMNormalClassicalMeasErrorDefinition("distanceToConspecific", 1, .1));
//		((MaximumLikelihoodEstimator) glmWithMeasError.getEstimator()).setLineSearchMethod(LineSearchMethod.HALF_STEP);
//		glmWithMeasError.doEstimation();
//		glmWithMeasError.getSummary();
//		Assert.assertTrue("Testing if convergence has been reached", glmWithMeasError.getEstimator().isConvergenceAchieved());
//		Assert.assertEquals("Testing the parameter estimate", -0.04971310978362395, glmWithMeasError.getParameters().getValueAt(1, 0), 1E-8);
//	}


 }
