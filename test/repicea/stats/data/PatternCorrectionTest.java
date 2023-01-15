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

import java.io.IOException;

import org.junit.Assert;
import org.junit.Ignore;
import org.junit.Test;

import repicea.io.javacsv.CSVReader;
import repicea.util.ObjectUtility;

public class PatternCorrectionTest {

	@Ignore
	@Test
	public void testStatusCorrection() throws IOException {
		String filename = ObjectUtility.getPackagePath(getClass()).replace("bin", "test") + "statusCorrected.csv";
		String refFilename = ObjectUtility.getPackagePath(getClass()).replace("bin", "test") + "statusCorrectedRef.csv";
		
		CSVReader reader = new CSVReader(filename);
		CSVReader refReader = new CSVReader(refFilename);
		Object[] record, refRecord;

		Assert.assertEquals("Testing number of records in each file", refReader.getRecordCount(), reader.getRecordCount()); 
		int i = 0;
		while ((record = reader.nextRecord()) != null) {
			refRecord = refReader.nextRecord();
			compareRecord(refRecord, record);
			i++;
		}
		reader.close();
		refReader.close();
		System.out.println("Number of records successfully tested (should be 1 676 351) = " + i);
	}

	private void compareRecord(Object[] refRecord, Object[] record) {
		Assert.assertEquals("Testing number of fields in each record", refRecord.length, record.length); 
		for (int i = 0; i < refRecord.length; i++) {
			Assert.assertEquals("Testing individual fields", refRecord[i], record[i]);
		}
	}
}
