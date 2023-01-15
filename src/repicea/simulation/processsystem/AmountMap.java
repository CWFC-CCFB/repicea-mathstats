/*
 * This file is part of the repicea-simulation library.
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
package repicea.simulation.processsystem;

import java.util.HashMap;
import java.util.Map;

/**
 * The AmountMap class contains all the quantities that can be contained in a ProcessUnit instance.
 * @author Mathieu Fortin - January 2014
 *
 * @param <E> an Enum class
 */
@SuppressWarnings("serial")
public class AmountMap<E extends Enum<?>> extends HashMap<E, Double> implements Cloneable {

	@Override
	public void putAll(Map<? extends E, ? extends Double> map) {
		for (E key : map.keySet()) {
			add(key, map.get(key));
		}
	}

	/**
	 * This method add a particular amount in a key.
	 * @param key
	 * @param amount
	 */
	public void add(E key, double amount) {
		if (containsKey(key)) {
			put(key, get(key) + amount);
		} else {
			put(key, amount);
		}
	}
	
	/**
	 * This method multiplies all the values in the current instance by a scalar.
	 * @param scalar the multiplicative factor
	 * @return a new AmountMap instance with updated values
	 */
	public AmountMap<E> multiplyByAScalar(double scalar) {
		AmountMap<E> outputMap = new AmountMap<E>();
		for (E key : keySet()) {
			outputMap.put(key, get(key) * scalar);
		}
		return outputMap;
	}

	@Override
	public AmountMap<E> clone() {
		AmountMap<E> outputMap = new AmountMap<E>();
		for (E key : keySet()) {
			outputMap.put(key, get(key));
		}
		return outputMap;
	}
	
	
	/**
	 * This method scales the AmountMap instances contained in more complex Map instances. It is recursive.
	 * IMPORTANT: the original map is changed. This avoid creating new instances which raise an exception
	 * when the Map class does not have any empty constructor.
	 * @param oMap
	 * @param scalar the multiplier factor
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void scaleMap(Map oMap, double scalar) {
//		try {
//			Map newMap = oMap.getClass().newInstance();
		for (Object key : oMap.keySet()) {
			Object value = oMap.get(key);
			if (value instanceof AmountMap) {
				oMap.put(key, ((AmountMap) value).multiplyByAScalar(scalar));
			} else {
				Map innerMap = (Map) value;
				scaleMap(innerMap, scalar);
			}
		}
//			return newMap;
//		} catch (Exception e) {
//			throw new InvalidParameterException("Unable to scale the Map instance!");
//		}
	}

//	/**
//	 * Adds the values to MonteCarloEstimate instances. If the map of MonteCarlo instances does not 
//	 * contain the key E, then the MonteCarloEstimate instance is created.
//	 * @param receivingMap
//	 */
//	public void addToMonteCarloEstimate(Map<E, MonteCarloEstimate> receivingMap) {
//		if (receivingMap == null) {
//			throw new InvalidParameterException("The receivingMap parameter cannot be null!");
//		}
//		Matrix value;
//		for (E element : keySet()) {
//			value = new Matrix(1,1);
//			value.m_afData[0][0] = get(element);
//			if (!receivingMap.containsKey(element)) {
//				receivingMap.put(element, new MonteCarloEstimate());
//			}
//			receivingMap.get(element).addRealization(value);
//		}
//	}

	
	
	
}
