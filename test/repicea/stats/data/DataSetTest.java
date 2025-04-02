/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2022 Mathieu Fortin for Rouge-Epicea.
 * Copyright (C) 2024 His Majesty the King in Right of Canada
 * Author: Mathieu Fortin, Canadian Forest Service
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

import java.io.IOException;
import java.lang.reflect.Array;
import java.math.BigInteger;
import java.util.Arrays;
import java.util.LinkedHashMap;

import org.junit.Assert;
import org.junit.Test;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;

public class DataSetTest {

	@Test
	public void test01CreateDataSetFromScratch() {
		DataSet myDataSet = new DataSet();
		myDataSet.addField("Field1", new Object[] {"true", "allo", "patate"});
		myDataSet.addField("Field2", new Object[] {"false", "hello", "carotte"});
		Object expected = "hello";
		Object observed = myDataSet.getObservations().get(1).values.get(1);
		Assert.assertEquals("Comparing values", observed, expected);
	}
	
	@SuppressWarnings("rawtypes")
	@Test
	public void test02JSONConversion() throws JsonProcessingException {
		DataSet myDataSet = new DataSet();
		myDataSet.addField("Field1", new Object[] {"true", "allo", "patate"});
		myDataSet.addField("Field2", new Object[] {"false", "hello", "carotte"});
		LinkedHashMap<String, Object>[] protoMap = myDataSet.getProtoMapArrayForJSONConversion();
		ObjectMapper mapper = new ObjectMapper();
		String jsonStr = mapper.writeValueAsString(protoMap);
		Object o = mapper.readValue(jsonStr, protoMap.getClass());
		Assert.assertTrue("Testing the instance type", o.getClass().isArray());
		Object entry1 = Array.get(o, 0);
		Assert.assertTrue("Testing first entry in the array", entry1 instanceof LinkedHashMap);
		String value1 = ((LinkedHashMap) entry1).get("Field1").toString();
		Assert.assertEquals("Testing value in Field1", "true", value1);
	}

	
	@Test
	public void test03LongInt() throws IOException {
		DataSet myDataSet = new DataSet(Arrays.asList(new String[] {"ID_PE", "TREENO"}));
		myDataSet.addObservation(new Object[]{"1000000000000", 4});
		Object idPeValue = myDataSet.getObservations().get(0).getValueAt(0);
		Assert.assertTrue("Testing BigInteger instance", idPeValue instanceof BigInteger);
		Assert.assertEquals("Testing string value of id pe", "1000000000000", idPeValue.toString());
	}

	@Test
	public void test04LongInt2() throws IOException {
		DataSet myDataSet = new DataSet(Arrays.asList(new String[] {"ID_PE", "TREENO"}));
		myDataSet.addObservation(new Object[]{"1000000000000", 4});
		myDataSet.addObservation(new Object[]{10, 3});
		myDataSet.indexFieldType();
		Class<?> clazz = myDataSet.getFieldTypes().get(0);
		Assert.assertTrue("Testing BigInteger instance", clazz.equals(BigInteger.class));
		for (Observation o : myDataSet.getObservations()) {
			Assert.assertTrue("Testing BigInteger instance", o.getValueAt(0) instanceof BigInteger);
		}
	}

	@Test
	public void test05DoubleParsing() {
		DataSet myDataSet = new DataSet();
		myDataSet.addField("Field1", new Object[] {"1.0", "1e-4", "2"});
		myDataSet.indexFieldType();
		Class<?> clazz = myDataSet.getFieldTypeOfThisField(myDataSet.getIndexOfThisField("Field1"));
		Assert.assertTrue("Testing Double instance", clazz.equals(Double.class));
	}

}
