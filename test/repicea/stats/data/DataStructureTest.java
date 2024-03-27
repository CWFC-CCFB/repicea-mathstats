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
package repicea.stats.data;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import org.junit.Assert;
import org.junit.BeforeClass;
import org.junit.Test;

import repicea.math.Matrix;
import repicea.stats.model.glm.GLModelTest;
import repicea.stats.model.lm.LinearModelTest;
import repicea.util.ObjectUtility;

public class DataStructureTest {

	private static DataSet DATASET;
	private static DataSet DATASET_ACC;
	
	@BeforeClass
	public static void initialization() throws Exception {
		DATASET = new DataSet(ObjectUtility.getPackagePath(GLModelTest.class) + "exampleDistanceSample.csv", true);
		DATASET_ACC = new DataSet(ObjectUtility.getPackagePath(LinearModelTest.class) + "diameterIncrementData.csv", true);
	}
	
	@Test
	public void simpleDataStructureTest() throws StatisticalDataException {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET);
		struct.setModelDefinition("isConspecificIn ~ distance");
		Matrix matrixX = struct.constructMatrixX();
		double intercept = matrixX.getValueAt(0, 0);
		Assert.assertEquals("Comparing intercept", 1d, intercept, 1E-8);
		double observedDistance = matrixX.getValueAt(0, 1);
		Assert.assertEquals("Comparing distances", 4.123105625617661, observedDistance, 1E-8);
	}

	@Test
	public void simpleDataStructureWithFormulaTest() throws StatisticalDataException {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET);
		struct.setModelDefinition("isConspecificIn ~ exp(distance)");
		Matrix matrixX = struct.constructMatrixX();
		double intercept = matrixX.getValueAt(0, 0);
		Assert.assertEquals("Comparing intercept", 1d, intercept, 1E-8);
		double observedDistance = matrixX.getValueAt(0, 1);
		Assert.assertEquals("Comparing distances", Math.exp(4.123105625617661), observedDistance, 1E-8);
	}

	@Test
	public void simpleDataStructureWithFormulaTest2() throws StatisticalDataException {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET);
		struct.setModelDefinition("isConspecificIn ~ exp(1 + distance) + exp(distance)");
		Matrix matrixX = struct.constructMatrixX();
		double intercept = matrixX.getValueAt(0, 0);
		Assert.assertEquals("Comparing intercept", 1d, intercept, 1E-8);
		double observedDistance = matrixX.getValueAt(0, 1);
		Assert.assertEquals("Comparing distances", Math.exp(1 + 4.123105625617661), observedDistance, 1E-8);
		double observedDistance2 = matrixX.getValueAt(0, 2);
		Assert.assertEquals("Comparing distances", Math.exp(4.123105625617661), observedDistance2, 1E-8);
	}
	
	@Test
	public void simpleDataStructureWithFormulaTest3() throws StatisticalDataException {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET);
		struct.setModelDefinition("isConspecificIn ~ exp(1 + distance) + log(distance)");
		Matrix matrixX = struct.constructMatrixX();
		double intercept = matrixX.getValueAt(0, 0);
		Assert.assertEquals("Comparing intercept", 1d, intercept, 1E-8);
		double observedDistance = matrixX.getValueAt(0, 1);
		Assert.assertEquals("Comparing distances", Math.exp(1 + 4.123105625617661), observedDistance, 1E-8);
		double observedDistance2 = matrixX.getValueAt(0, 2);
		Assert.assertEquals("Comparing distances", Math.log(4.123105625617661), observedDistance2, 1E-8);
	}

	@Test
	public void simpleDataStructureWithFormulaTest4() throws StatisticalDataException {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET);
		try {
			struct.setModelDefinition("isConspecificIn ~ exp(1 + distance) + log()");
			Assert.fail("The setModelDefinition method was supposed to throw an Exception!");
		} catch (InvalidParameterException e) {}
	}

	@Test
	public void hierarchicalDataStructureWithFormulaTest1() throws StatisticalDataException {
		HierarchicalStatisticalDataStructure struct = new GenericHierarchicalStatisticalDataStructure(DATASET);
		struct.setModelDefinition("isConspecificIn ~ exp(1 + distance) + exp(distance) + (distance | id)");
		Matrix matrixX = struct.constructMatrixX();
		double intercept = matrixX.getValueAt(0, 0);
		Assert.assertEquals("Comparing intercept", 1d, intercept, 1E-8);
		double observedDistance = matrixX.getValueAt(0, 1);
		Assert.assertEquals("Comparing distances", Math.exp(1 + 209.69024774652732), observedDistance, 1E-8);
		double observedDistance2 = matrixX.getValueAt(0, 2);
		Assert.assertEquals("Comparing distances", Math.exp(209.69024774652732), observedDistance2, 1E-8);
	}

	@Test
	public void inclusiveInteractionDecompositionTwoEffects() {
		List<String> twoEffects = new ArrayList<String>();
		twoEffects.add("patate");
		twoEffects.add("carotte");
		List<String> allComb = GenericStatisticalDataStructure.getAllCombinations(twoEffects);
		Assert.assertEquals("Testing number of interactions", 1, allComb.size());
	}

	@Test
	public void inclusiveInteractionDecompositionThreeEffects() {
		List<String> threeEffects = new ArrayList<String>();
		threeEffects.add("patate");
		threeEffects.add("carotte");
		threeEffects.add("navet");
		
		List<String> allComb = GenericStatisticalDataStructure.getAllCombinations(threeEffects);
		Assert.assertEquals("Testing number of interactions", 4, allComb.size());
	}

	
	@Test
	public void genericDataStructureWithFormulaTest1() throws Exception {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET_ACC);
		struct.setModelDefinition("yTrans ~ BAL*dbhCm.x");
		Assert.assertTrue("Testing number of effects", struct.effects.size() == 3);
		Assert.assertTrue("Testing BAL", struct.effects.containsKey("BAL"));
		Assert.assertTrue("Testing dbhCm.x", struct.effects.containsKey("dbhCm.x"));
		Assert.assertTrue("Testing BAL:dbhCm.x", struct.effects.containsKey("BAL:dbhCm.x"));
	}

	@Test
	public void genericDataStructureWithFormulaTest2() throws Exception {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET_ACC);
		struct.setModelDefinition("yTrans ~ BAL*dbhCm.x + dbhCm.x");
		Assert.assertTrue("Testing number of effects", struct.effects.size() == 3);
		Assert.assertTrue("Testing BAL", struct.effects.containsKey("BAL"));
		Assert.assertTrue("Testing dbhCm.x", struct.effects.containsKey("dbhCm.x"));
		Assert.assertTrue("Testing BAL:dbhCm.x", struct.effects.containsKey("BAL:dbhCm.x"));
	}

	@Test
	public void genericDataStructureWithFormulaTest3() throws Exception {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET_ACC);
		struct.setModelDefinition("yTrans ~ BAL*dbhCm.x + dbhCm.x:BAL");
		Assert.assertTrue("Testing number of effects", struct.effects.size() == 3);
		Assert.assertTrue("Testing BAL", struct.effects.containsKey("BAL"));
		Assert.assertTrue("Testing dbhCm.x", struct.effects.containsKey("dbhCm.x"));
		Assert.assertTrue("Testing BAL:dbhCm.x", struct.effects.containsKey("BAL:dbhCm.x"));
	}

	@Test
	public void genericDataStructureWithFormulaTest4() throws Exception {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET_ACC);
		struct.setModelDefinition("yTrans ~ dbhCm.x*BAL");
		Assert.assertTrue("Testing number of effects", struct.effects.size() == 3);
		Assert.assertTrue("Testing BAL", struct.effects.containsKey("BAL"));
		Assert.assertTrue("Testing dbhCm.x", struct.effects.containsKey("dbhCm.x"));
		Assert.assertTrue("Testing BAL:dbhCm.x", struct.effects.containsKey("BAL:dbhCm.x"));
	}

	@Test
	public void genericDataStructureWithFormulaTest5() throws Exception {
		GenericStatisticalDataStructure struct = new GenericStatisticalDataStructure(DATASET_ACC);
		struct.setModelDefinition("yTrans ~ exp(BAL)*dbhCm.x");
		Assert.assertTrue("Testing number of effects", struct.effects.size() == 3);
		Assert.assertTrue("Testing BAL", struct.effects.containsKey("&0"));
		Assert.assertTrue("Testing dbhCm.x", struct.effects.containsKey("dbhCm.x"));
		Assert.assertTrue("Testing BAL:dbhCm.x", struct.effects.containsKey("&0:dbhCm.x"));
		Matrix matX = struct.constructMatrixX();
		Assert.assertEquals("Nb cols in X", 4, matX.m_iCols);
		double BAL = (Double) DATASET_ACC.getValueAt(0, "BAL");
		double dbhCm = (Double) DATASET_ACC.getValueAt(0, "dbhCm.x");
		Assert.assertEquals("Value of intercept", 1d, matX.getValueAt(0, 0), 1E-8);
		Assert.assertEquals("Value of exp(BAL)", Math.exp(BAL), matX.getValueAt(0, 1), 1E-8);
		Assert.assertEquals("Value of dbhCm.x", dbhCm, matX.getValueAt(0, 2), 1E-8);
		Assert.assertEquals("Value of exp(BAL):dbhCm.x", Math.exp(BAL)*dbhCm, matX.getValueAt(0, 3), 1E-8);
	}

	// TODO add a test that has an interaction between an effect and a formula
}
