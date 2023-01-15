/*
 * This file is part of the repicea-simulation library.
 *
 * Copyright (C) 2009-2017 Mathieu Fortin for Rouge-Epicea
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

import java.util.ArrayList;
import java.util.List;

@SuppressWarnings("rawtypes")
public class TestProcessUnit extends ProcessUnit {

	protected final List<Processor> processorList;
	boolean partOfEndlessLoop;
	
	protected TestProcessUnit() {
		processorList = new ArrayList<Processor>();
	}

//	protected TestProcessUnit(List<Processor> processorList) {
//		this();
//		this.processorList.addAll(processorList);
//	}

	protected TestProcessUnit createNewProcessUnitFromThisOne() {
		TestProcessUnit tpu = new TestProcessUnit();
		tpu.processorList.addAll(processorList);
		return tpu;
	}
	
	protected boolean recordProcessor(Processor processor) {
		partOfEndlessLoop = processorList.contains(processor);
		processorList.add(processor);
		if (partOfEndlessLoop) {
			int indexBegin = processorList.indexOf(processorList.get(processorList.size() - 1));
			for (int i = indexBegin; i < processorList.size(); i++) {
				processorList.get(i).setPartOfEndlessLoop(true);
			}
		}
		return partOfEndlessLoop;
	}
	
//	protected void reset() {
//		partOfEndlessLoop = false;
//		processorList.clear();
//	}
}
