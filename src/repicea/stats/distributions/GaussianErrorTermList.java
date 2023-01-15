/*
 * This file is part of the repicea-statistics library.
 *
 * Copyright (C) 2009-2014 Mathieu Fortin for Rouge-Epicea
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
package repicea.stats.distributions;

import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import repicea.math.Matrix;

@SuppressWarnings("serial")
public final class GaussianErrorTermList extends ArrayList<GaussianErrorTerm> {

	protected boolean updated;
	
	/**
	 * This interface ensures the instance can return an index that will serve as distance for
	 * the calculation of the variance-covariance matrix.
	 * @author Mathieu Fortin - August 2014
	 */
	public static interface IndexableErrorTerm {
		
		/**
		 * This method returns the index of the error term. Typically, this value is the time and it serves to 
		 * calculate the distance between two observations when computing the variance-covariance matrix.
		 * @return an Integer
		 */
		public int getErrorTermIndex();
	}
	
	public List<Integer> getDistanceIndex() {
		List<Integer> indexList = new ArrayList<Integer>();
		for (GaussianErrorTerm res : this) {
			indexList.add(res.distanceIndex);
		}
		return indexList;
	}

	public Matrix getNormalizedErrors() {
		Matrix mat = new Matrix(size(),1);
		for (int i = 0; i < size(); i++) {
			mat.setValueAt(i, 0, get(i).normalizedValue);
		}
		return mat;
	}

	protected Matrix getRealizedErrors() {
		Matrix mat = new Matrix(size(),1);
		for (int i = 0; i < size(); i++) {
			mat.setValueAt(i, 0, get(i).value);
		}
		return mat;
	}

	public void updateErrorTerm(Matrix errorTerms) {
		for (int i = 0; i < errorTerms.m_iRows; i++) {
			GaussianErrorTerm error = get(i);
			if (error.value == null) {
				error.value = errorTerms.getValueAt(i, 0);
			}
		}
		updated = true;
	}
	
	public double getErrorForIndexableInstance(IndexableErrorTerm indexableErrorTerm) {
		int distanceIndex = indexableErrorTerm.getErrorTermIndex();
		int index = getDistanceIndex().indexOf(distanceIndex);
		if (index < 0) {
			throw new InvalidParameterException("This distance index is not contained in the GaussianErrorTermList");
		} else {
			return get(index).value;
		}
	}

	@Override
	public boolean add(GaussianErrorTerm term) {
		boolean result = super.add(term);
		updated = false;
		Collections.sort(this);
		return result;
	}
	
	
	
//	protected static class FakeClass implements IndexableErrorTerm {
//
//		int index;
//		
//		FakeClass(int index) {
//			this.index = index;
//		}
//		
//		@Override
//		public int getErrorTermIndex() {
//			return index;
//		}
//		
//	}
//	
//	
//	public static void main(String[] args) {
//		GaussianErrorTermList list = new GaussianErrorTermList();
//		list.add(new GaussianErrorTerm(new FakeClass(2009)));
//		list.add(new GaussianErrorTerm(new FakeClass(2029)));
//		list.add(new GaussianErrorTerm(new FakeClass(2019)));
//	}
}

