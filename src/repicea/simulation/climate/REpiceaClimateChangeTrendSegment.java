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
package repicea.simulation.climate;

import java.io.Serializable;

class REpiceaClimateChangeTrendSegment implements Serializable {

	final int startDateYr;
	final int endDateYr;
	final REpiceaClimateVariableChangeMap changeMap;
	
	REpiceaClimateChangeTrendSegment(int startDateYr, int endDateYr, REpiceaClimateVariableChangeMap changeMap) {
		this.startDateYr = startDateYr;
		this.endDateYr = endDateYr;
		this.changeMap = new REpiceaClimateVariableChangeMap();
		this.changeMap.putAll(changeMap); // the change map is cloned to avoid future concurrent modification
	}
	
}
