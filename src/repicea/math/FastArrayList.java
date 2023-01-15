/*
 * This file is part of the repicea-statistics library.
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
package repicea.math;

import java.io.Serializable;

/**
 * The FastArrayList class is used by any AbstractMathematicalFunction-derived class. It is instantiated by and 
 * only by this class. If the List contains Number instances, the getMatrix() function can return the list into a
 * Matrix instance.
 * @author Mathieu Fortin - November 2012
 * @param <T> the class of the object in this list
 */
@SuppressWarnings("serial")
public class FastArrayList<T extends Serializable> implements Serializable {

	private int capacity = 10;
	private int length = 0;
	
	private Object[] internalArray;
	
	FastArrayList() {
		internalArray = new Object[10];
	}

	public void add(T obj) {
		if (length == capacity) {
			Object[] copyInternal = internalArray;
			capacity += 10;
			internalArray = new Object[capacity];
			for (int i = 0; i < length; i++) {
				internalArray[i] = copyInternal[i];
			}
		}
		internalArray[length] = obj;
		length++;
	}

	public boolean contains(T obj) {
		return indexOf(obj) > -1;
	}

	public void clear() {
		length = 0;
	}

	public int indexOf(T obj) {
		for (int i = 0; i < length; i++) {
			if (internalArray[i].equals(obj)) {
				return i;
			}
		}
		return -1;
	}
	
	public int size() {return length;}
	
	public void set(int i, T obj) {
		internalArray[i] = obj;
	}
	
	@SuppressWarnings("unchecked")
	public T get(int i) {
		if (i >= length) {
			return null;
		} else {
			return (T) internalArray[i];
		}
	}
	
	public Matrix getMatrix() {
		if (length > 0 && internalArray[0] instanceof Number) {
			Matrix outputMatrix = new Matrix(length, 1);
			for (int i = 0; i < length; i++) {
				outputMatrix.setValueAt(i, 0, ((Number) internalArray[i]).doubleValue());
			}
			return outputMatrix;
		} else {
			return null;
		}
	} 
	
	@SuppressWarnings("unchecked")
	protected T[] getInternalArray() {return (T[]) internalArray;}
	
}
