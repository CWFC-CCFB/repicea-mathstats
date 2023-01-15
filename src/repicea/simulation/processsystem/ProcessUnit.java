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

import java.util.Map;


/**
 * The ProcessUnit class represents a basic unit in a process system.
 * @author Mathieu Fortin - January 2014
 *
 * @param <E> an enum instance
 */
public class ProcessUnit<E extends Enum<?>> {

	private AmountMap<E> amountMap;
	
	protected ProcessUnit() {
		amountMap = new AmountMap<E>();
	}
	
	protected ProcessUnit(Map<E, Double> amountMap) {
		this();
		this.amountMap.putAll(amountMap);
	}
	
	/**
	 * This method returns the quantities contained in this process unit.
	 * @return an AmountMap instance
	 */
	public AmountMap<E> getAmountMap() {
		return amountMap;
	}
	
	protected void addProcessUnit(ProcessUnit<E> unit) {
		amountMap.putAll(unit.amountMap);
	}

	
	
	
}
