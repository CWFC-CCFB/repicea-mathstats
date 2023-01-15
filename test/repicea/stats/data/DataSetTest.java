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

import org.junit.Assert;
import org.junit.Test;


public class DataSetTest {

	@Test
	public void createDataSetFromScratch() {
		DataSet myDataSet = new DataSet();
		myDataSet.addField("Field1", new Object[] {"true", "allo", "patate"});
		myDataSet.addField("Field2", new Object[] {"false", "hello", "carotte"});
		Object expected = "hello";
		Object observed = myDataSet.getObservations().get(1).values.get(1);
		Assert.assertEquals("Comparing values", observed, expected);
	}
}
