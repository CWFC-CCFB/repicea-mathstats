/*
 * This file is part of the repicea library.
 *
 * Copyright (C) 2009-2018 Mathieu Fortin for Rouge-Epicea
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
package repicea.simulation.covariateproviders.plotlevel;

/**
 * This interface ensures that a particular plot can provide its composition, that is
 * whether it is broadleaved dominated, conifer dominated or mixed.
 * @author Mathieu Fortin - May 2018
 *
 */
public interface SpeciesCompositionProvider {

	public static enum SpeciesComposition {
		BroadleavedDominated,
		ConiferDominated,
		Mixed
	}
	
	/**
	 * This method returns the composition of the plot, that is either broadleaved dominated, 
	 * conifer dominated or mixed.
	 * @return a SpeciesComposition enum
	 */
	public SpeciesComposition	getSpeciesComposition();
	
}
