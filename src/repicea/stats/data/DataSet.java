/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2019 Mathieu Fortin for Rouge-Epicea
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
import java.security.InvalidParameterException;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import com.fasterxml.jackson.annotation.JsonIgnore;

import repicea.data.Record;
import repicea.data.Table;
import repicea.gui.REpiceaUIObject;
import repicea.gui.components.REpiceaTable;
import repicea.gui.components.REpiceaTableModel;
import repicea.io.FormatField;
import repicea.io.FormatReader;
import repicea.io.FormatWriter;
import repicea.io.GExportFieldDetails;
import repicea.io.Saveable;
import repicea.math.Matrix;

/**
 * The DataSet class contains many observations and implements the method to read a dataset with a FormatReader instance.
 * @author Mathieu Fortin - November 2012
 */
public class DataSet implements Table, Saveable, REpiceaUIObject {

//	private static enum MessageID implements TextableEnum {
//
//		ReadingFileMessage("Reading file...", "Lecture du fichier...");
//
//		MessageID(String englishText, String frenchText) {
//			setText(englishText, frenchText);
//		}
//		
//		@Override
//		public void setText(String englishText, String frenchText) {
//			REpiceaTranslator.setString(this, englishText, frenchText);
//		}
//
//		@Override
//		public String toString() {return REpiceaTranslator.getString(this);}
//	}
	

	
	protected List<String> fieldNames;
	protected List<Class<?>> fieldTypes;
	protected List<Observation> observations;
	private final String originalFilename;
	
	@JsonIgnore
	private transient REpiceaTable table;
	@JsonIgnore
	private Map<Integer, NumberFormat> formatters;
	
	protected DataSet(String filename) {
		this.originalFilename = filename;
		fieldNames = new ArrayList<String>();
		fieldTypes = new ArrayList<Class<?>>();
		observations = new ArrayList<Observation>();
	}

	/**
	 * General constructor.
	 * @param filename the name of the file to be read
	 * @param autoLoad true if the file is to be read now. Typically, this boolean is set to false when the swingworker is
	 * launched from a window that retrieves some events.
	 * @throws Exception if an error has occurred during the loading
	 */
	public DataSet(String filename, boolean autoLoad) throws Exception {
		this(filename);
		if (autoLoad) {
			load();
		}
	}
	
	/**
	 * An empty DataSet instance with known field names.
	 * @param fieldNames a List of String instances that represent the field names
	 */
	public DataSet(List<String> fieldNames) {
		this();
		for (String fieldName : fieldNames) {
			addFieldName(fieldName);
		}
	}

	/**
	 * An empty DataSet to be populated using the addField method.
	 */
	public DataSet() {
		this((String) null);
	}

	/**
	 * This method returns any object in the dataset at row i and column j.
	 * @param i the index of the row
	 * @param j the index of the column
	 * @return an Object instance
	 */
	protected Object getValueAt(int i, int j) {
		return observations.get(i).values.get(j);
	}

	public Object getValueAt(int i, String fieldName) {
		int j = getIndexOfThisField(fieldName);
		if (j != -1) {
			return getValueAt(i,j);
		} else {
			return null;
		}
	}

	protected final void setValueAt(int i, int j, Object value) {
		if (value.getClass().equals(fieldTypes.get(j))) {
			observations.get(i).values.remove(j);
			observations.get(i).values.add(j, value);
		} 
	}
	
	/**
	 * Indexes the different field types. More specifically, it goes 
	 * through the columns and find the appropriate class for a particular
	 * field. This method should be called after adding all the observations.
	 */
	public void indexFieldType() {
		fieldTypes.clear();
		for (int j = 0; j < fieldNames.size(); j++) {
			setClassOfThisField(j);
		}
	}
	
	private void setClassOfThisField(int fieldIndex) {
		if (isInteger(fieldIndex)) {
			setFieldType(fieldIndex, Integer.class);
		} else if (isDouble(fieldIndex)) {
			setFieldType(fieldIndex, Double.class);
			reconvertToDoubleIfNeedsBe(fieldIndex);
		} else {
			setFieldType(fieldIndex, String.class);
			reconvertToStringIfNeedsBe(fieldIndex);
		}
	}

	private void setFieldType(int fieldIndex, Class<?> clazz) {
		if (fieldIndex < fieldTypes.size()) {
			fieldTypes.set(fieldIndex, clazz);	
		} else if (fieldIndex == fieldTypes.size()) {
			fieldTypes.add(clazz);	
		} else {
			throw new InvalidParameterException("The field type cannot be set!");
		}
	}
		
	private boolean isInteger(int j) {
		boolean isInteger = true;
		for (int i = 0; i < getNumberOfObservations(); i++) {
			if (!(getValueAt(i,j) instanceof Integer)) {
					isInteger = false;
					break;
			} 
		}
		return isInteger;
	}


	private void reconvertToStringIfNeedsBe(int j) {
		for (int i = 0; i < getNumberOfObservations(); i++) {
			Object value = getValueAt(i,j);
			if ((value instanceof Number)) {
				setValueAt(i,j, value.toString());
			}
		} 
	}

	private void reconvertToDoubleIfNeedsBe(int j) {
		for (int i = 0; i < getNumberOfObservations(); i++) {
			Object value = getValueAt(i,j);
			if ((value instanceof Integer)) {
				setValueAt(i,j, ((Integer) value).doubleValue()); // MF2020-04-30 Bug corrected 
			}
		} 
	}


	private boolean isDouble(int indexJ) {
		boolean isDouble = true;
		for (int i = 0; i < getNumberOfObservations(); i++) {
			if (!(getValueAt(i,indexJ) instanceof Number)) {
					isDouble = false;
					break;
			} 
		}
		return isDouble;
	}
	
	/**
	 * This method returns the index of a particular field.
	 * @param fieldName the name of the field
	 * @return an integer
	 */
	public int getIndexOfThisField(String fieldName) {return fieldNames.indexOf(fieldName);}

	/**
	 * This method sorts the data according to the fields represented by the indices in fieldIndices parameter.
	 * @param fieldIndices a List of field indices
	 */
	@SuppressWarnings("unchecked")
	public void sortObservations(List<Integer> fieldIndices) {
		Observation.comparableFields = fieldIndices;
		Collections.sort(observations);
	}
	
	
	/**
	 * This method returns the number of observations in the dataset.
	 * @return an integer
	 */
	public int getNumberOfObservations() {
		return observations.size();
	}
	
	
	@SuppressWarnings("rawtypes")
	protected Class getFieldTypeOfThisField(int i) {
		return fieldTypes.get(i);
	}
	
	@SuppressWarnings("rawtypes")
	protected Class getFieldTypeOfThisField(String fieldName) {
		int index = getIndexOfThisField(fieldName);
		if (index == -1) 
			throw new InvalidParameterException("Field " + fieldName + " cannot be found in the dataset!");
		else 
			return getFieldTypeOfThisField(getIndexOfThisField(fieldName));
	}
	

	protected Matrix getVectorOfThisField(String fieldName) {
		return getVectorOfThisField(getIndexOfThisField(fieldName));
	}


	protected Matrix getVectorOfThisField(int j) {
		Matrix output = new Matrix(observations.size(), 1);
		for (int i = 0; i < observations.size(); i++) {
			output.setValueAt(i, 0, ((Number) getValueAt(i,j)).doubleValue());
		}
		return output;
	}
	
	
	
	
	/**
	 * This method returns all the possible values in this field.  
	 * @param j the index of the field.
	 * @return a SORTED list of all the possible value.
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	protected List getPossibleValuesInThisField(int j) {
		List possibleValues = new ArrayList();
		for (int i = 0; i < observations.size(); i++) {
			Object value = getValueAt(i,j);
			if (!possibleValues.contains(value)) {
				possibleValues.add(value);
			}
		}
		Collections.sort(possibleValues);
		return possibleValues;
	}

	protected Matrix getDummyMatrix(List<?> possibleValues, int fieldIndex) {
		Matrix outputMatrix = new Matrix(getNumberOfObservations(), possibleValues.size());
		for (int i = 0; i < getNumberOfObservations(); i++) {
			int position = possibleValues.indexOf(getValueAt(i, fieldIndex));
			if (position >= 0 && position < outputMatrix.m_iCols) {
				outputMatrix.setValueAt(i, position, 1d);
			}
		}
		return outputMatrix;
	}
	
	public void addObservation(Object[] observationFrame) {
		parseDifferentFields(observationFrame);
		observations.add(new Observation(observationFrame));
	}
	
	private void addFieldName(String originalName) {
		int index = 0;
		String name = originalName;
		while (fieldNames.contains(name)) {
			name = originalName + index;
			index++;
		}
		fieldNames.add(name);
	}
	
	public void addField(String name, Object[] field) {
		if (observations.size() > 0 && field.length != observations.size()) {	// will only trigger if there are some observations already
			throw new InvalidParameterException("The number of observations in the new field does not match the number of observations in the dataset!");
		}
		addFieldName(name);
		
		for (int i = 0; i < field.length; i++) {
			if (i < observations.size()) { // means the observation exists already 
				observations.get(i).values.add(field[i]);
			} else {
				observations.add(new Observation(new Object[] {field[i]}));
			}
		}
		
		setClassOfThisField(fieldNames.size() - 1);
//		if (isDouble(fieldNames.size() - 1)) {
//			fieldTypes.add(Double.class);
//		} else {
//			fieldTypes.add(String.class);
//		}
	}

	@Override
	public void save(String filename) throws IOException {
		try {
			FormatWriter<?> writer = FormatWriter.createFormatWriter(false, filename);
			GExportFieldDetails exportField;
			List<FormatField> headerFields = new ArrayList<FormatField>();
			Object[] record;
			for (int i = 0; i < observations.size(); i++) {
				record = new Object[fieldNames.size()];
				for (int j = 0; j < fieldNames.size(); j++) {
					record[j] = getValueAt(i,j);
					if (i == 0) {
						exportField = new GExportFieldDetails(fieldNames.get(j), getValueAt(i,j));
						headerFields.add(writer.convertGExportFieldDetailsToFormatField(exportField));
					}
				}
				if (i == 0) {
					writer.setFields(headerFields);
				}
				writer.addRecord(record);
			}
			writer.close();
		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	
	private void load() throws Exception {
		fieldNames.clear();
		observations.clear();

		try {
			FormatReader<?> reader = FormatReader.createFormatReader(originalFilename);
			FormatField field;
			for (int i = 0; i < reader.getHeader().getNumberOfFields(); i++) {
				field = reader.getHeader().getField(i);
				fieldNames.add(field.getName());
			}

//			int nbRecords = reader.getRecordCount();
//			int recordsRead = 0;
			
//			firePropertyChange(REpiceaProgressBarDialog.LABEL, 0d, MessageID.ReadingFileMessage.toString());
			
			Object[] lineRead = reader.nextRecord();
			while (lineRead != null) {
				addObservation(lineRead);
//				recordsRead++;
//				int progress = (int) ((recordsRead * 100d) / nbRecords);
//				firePropertyChange(REpiceaProgressBarDialog.PROGRESS, recordsRead, progress);
				lineRead = reader.nextRecord();
			}
			
			indexFieldType();
			

		} catch (IOException e1) {
			e1.printStackTrace();
		}
	}

	private void parseDifferentFields(Object[] lineRead) {
		for (int i = 0; i < fieldNames.size(); i++) {
			try {
				lineRead[i] = (Integer) Integer.parseInt(lineRead[i].toString());
			} catch (NumberFormatException e1) {
				try {
					lineRead[i] = (Double) Double.parseDouble(lineRead[i].toString());
				} catch (NumberFormatException e2) {
					lineRead[i] = lineRead[i].toString();
				}
			}
		}
	}
		
	@Override
	public List<String> getFieldNames() {
		List<String> fieldNames = new ArrayList<String>();
		fieldNames.addAll(this.fieldNames);
		return fieldNames;
	}
	
	@Override
	public List<Class<?>> getFieldTypes() {
		List<Class<?>> fieldTypes = new ArrayList<Class<?>>();
		fieldTypes.addAll(this.fieldTypes);
		return fieldTypes;
	}
	
	@JsonIgnore
	@Override
	public REpiceaTable getUI() {
		if (table == null) {
			REpiceaTableModel model = new REpiceaTableModel(this);
			table = new REpiceaTable(model, false);	// no pop up 
		}
		return table;
	}


	@JsonIgnore
	@Override
	public boolean isVisible() {
		return table.isVisible();
	}


	/**
	 * Returns the observations of the data set.
	 * @return a List of Observation instances
	 */	
	public List<Observation> getObservations() {
		return observations;
	}
	
	
	private static String addSpacesUpTo(String str, int desiredLength, int buffer, boolean isNumber) {
		StringBuilder sb = new StringBuilder();
		if (isNumber) {
			int nbSpaceToAdd = desiredLength - str.length();
			while (sb.length() < nbSpaceToAdd)
				sb.append(" ");
			sb.append(str);
		} else {
			sb.append(str);
			while (sb.length() < desiredLength)
				sb.append(" ");
		}
		while(sb.length() < desiredLength + buffer) 
			sb.append(" ");
		return sb.toString();
	}
	
	@JsonIgnore
	private Map<Integer, NumberFormat> getFormatters() {
		if (formatters == null) {
			formatters = new HashMap<Integer, NumberFormat>();
		}
		return formatters;
	}
	
	@JsonIgnore
	private NumberFormat getFormatter(int j) {
		return getFormatters().get(j);
	}
	
	/**
	 * Set a NumberFormat instance for a particular field
	 * @param fieldId the id of the field
	 * @param formatter a NumberFormat instance
	 */
	public void setFormatter(int fieldId, NumberFormat formatter) {
		if (fieldId < 0 || fieldId >= this.getFieldNames().size()) {
			throw new InvalidParameterException("The fieldId argument must be positive (>= 0) and smaller than the number of fields!");
		}
		getFormatters().put(fieldId, formatter);
	}

	/**
	 * Clear the NumberFormat instances associated with the fields.
	 */
	public void clearFormatters() {
		getFormatters().clear();
	}
	
	@Override
	public String toString() {
		boolean exceeds = false;
		int maxLength;
		if (getNumberOfObservations() <= 100) {
			maxLength = getNumberOfObservations();
		} else {
			exceeds = true;
			maxLength = 100;
		}
		int nbFields = getFieldNames().size();
		int[] maxSizes = new int[nbFields];
		for (int j = 0; j < nbFields; j++)  {
			String fieldName = getFieldNames().get(j);
			maxSizes[j] = fieldName.length();
		}
		
		for (int i = 0; i < maxLength; i++) {
			for (int j = 0; j < nbFields; j++)  {
				Object o = getObservations().get(i).getValueAt(j);
				int currentLength = o instanceof Double && getFormatter(j) != null ? getFormatter(j).format((Double) o).length() : o.toString().length();
				if (maxSizes[j] < currentLength) {
					maxSizes[j] = currentLength;
				}
			}
		}
		StringBuilder output = new StringBuilder();
		for (int j = 0; j < nbFields; j++)  {
			String fieldName = getFieldNames().get(j);
			output.append(addSpacesUpTo(fieldName, maxSizes[j], 2, false));	// not number instances
		}
		
		output.append(System.lineSeparator());
		
		for (int i = 0; i < maxLength; i++) {
			for (int j = 0; j < nbFields; j++)  {
				Object o = getObservations().get(i).getValueAt(j);
				String value = o instanceof Double && getFormatter(j) != null ? getFormatter(j).format((Double) o) : o.toString();
				output.append(addSpacesUpTo(value, maxSizes[j], 2, o instanceof Number));
			}
			output.append(System.lineSeparator());
		}
		if (exceeds) {
			output.append("Only " + maxLength + " out of " + getNumberOfObservations() + " observations printed!");
		}
		return output.toString();
	}

	/**
	* This method returns a list of the values in a particular field.
	* @param i the field id
	* @return a List of object instance
	*/	
	public List<Object> getFieldValues(int i) {
		List<Object> objs = new ArrayList<Object>();
		for (Observation obs : observations) {
			objs.add(obs.values.get(i));
		}
		return objs;
	}
	
	@Override
	public List<? extends Record> getRecords() {
		return getObservations();
	}
}
