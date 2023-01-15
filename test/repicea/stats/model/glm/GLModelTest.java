/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2012 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.model.glm;

import static org.junit.Assert.assertEquals;

import java.security.InvalidParameterException;
import java.util.logging.ConsoleHandler;
import java.util.logging.Level;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import repicea.math.optimizer.NewtonRaphsonOptimizer;
import repicea.stats.data.DataSet;
import repicea.stats.estimators.MaximumLikelihoodEstimator;
import repicea.stats.model.glm.Family.GLMDistribution;
import repicea.stats.model.glm.LinkFunction.Type;
import repicea.stats.model.glm.copula.FGMCopulaGLModelTest;
import repicea.util.ObjectUtility;
import repicea.util.REpiceaLogManager;

public class GLModelTest {

	@BeforeClass
	public static void doThisBefore() {
		Level l = Level.WARNING;
		NewtonRaphsonOptimizer.LOGGER_NAME = MaximumLikelihoodEstimator.LOGGER_NAME;
		ConsoleHandler ch = new ConsoleHandler();
		ch.setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).setLevel(l);
		REpiceaLogManager.getLogger(MaximumLikelihoodEstimator.LOGGER_NAME).addHandler(ch);		
	}

	@Test
    public void TestWithSimpleGLModel() throws Exception {
		double expectedLlk = -1091.9193286646055;
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.Logit, "coupe ~ diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
		glm.doEstimation();
		double actualLlk = glm.getCompleteLogLikelihood().getValue();
		assertEquals(expectedLlk, actualLlk, 1E-5);
	}

	@Test
    public void TestWithOffset() throws Exception {
		String filename = ObjectUtility.getPackagePath(FGMCopulaGLModelTest.class).concat("donneesR_min.csv");
		DataSet dataSet = new DataSet(filename, true);
		try {
			new GeneralizedLinearModel(dataSet, GLMDistribution.Bernoulli, Type.Logit, "coupe ~ offset(TIGES_ID) + diffdhp + marchand:diffdhp + marchand:diffdhp2 +  essence");
			Assert.fail("The GeneralizedLinearModel instance should have thrown an Exception!");
		} catch(InvalidParameterException e) {}
	}

	@Test
    public void TestWithNegBinomial() throws Exception {
		String filename = ObjectUtility.getPackagePath(GLModelTest.class).concat("recruitERS.csv");
		DataSet dataSet = new DataSet(filename, true);
		GeneralizedLinearModel glm = new GeneralizedLinearModel(dataSet, GLMDistribution.NegativeBinomial, Type.Log, "y ~ dt + G_F + G_R + speciesThere");
		glm.doEstimation();
		System.out.println(glm.getSummary());
		Assert.assertTrue("Testing if convergence was achieved", glm.getEstimator().isConvergenceAchieved());
		double parmEst_intercept = glm.getParameters().getValueAt(0, 0);
		double expected = 0.250792010984403;
		Assert.assertEquals("Testing parameter estimate: intercept", expected, parmEst_intercept, 1E-8);
		double parmEst_dt = glm.getParameters().getValueAt(1, 0);
		expected = 0.06509727794589154;
		Assert.assertEquals("Testing parameter estimate: dt", expected, parmEst_dt, 1E-8);
		double parmEst_GF = glm.getParameters().getValueAt(2, 0);
		expected = -0.05175664359330031;
		Assert.assertEquals("Testing parameter estimate: G_F", expected, parmEst_GF, 1E-8);
		double parmEst_GR = glm.getParameters().getValueAt(3, 0);
		expected = -0.11805966012849123;
		Assert.assertEquals("Testing parameter estimate: G_R", expected, parmEst_GR, 1E-8);
		double parmEst_speciesThere = glm.getParameters().getValueAt(4, 0);
		expected = 0.4911558058634061;
		Assert.assertEquals("Testing parameter estimate: speciesThere", expected, parmEst_speciesThere, 1E-8);
		double parmEst_dispersion = glm.getParameters().getValueAt(5, 0);
		expected = 1.2474827482825246;
		Assert.assertEquals("Testing parameter estimate: dispersion", expected, parmEst_dispersion, 1E-8);
		
	}

}
