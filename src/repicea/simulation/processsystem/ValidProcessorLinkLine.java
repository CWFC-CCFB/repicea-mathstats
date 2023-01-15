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

import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ComponentEvent;
import java.awt.event.ComponentListener;

@SuppressWarnings("serial")
public abstract class ValidProcessorLinkLine extends AbstractProcessorLinkLine implements ComponentListener {

	
	protected ValidProcessorLinkLine(SystemPanel panel,	AnchorProvider fatherAnchor, AnchorProvider sonAnchor) {
		super(panel, fatherAnchor, sonAnchor);
		setBackground(Color.BLACK);
		Dimension dim = new Dimension(8,8);
		setMinimumSize(dim);
		setPreferredSize(dim);
		setMaximumSize(dim);
		panel.add(this);
	}

	@Override
	public void componentHidden(ComponentEvent arg0) {}

	@Override
	public void componentMoved(ComponentEvent arg0) {
		panel.repaint();
	}

	@Override
	public void componentResized(ComponentEvent arg0) {}
	
	@Override
	public void componentShown(ComponentEvent arg0) {}


	@Override
	protected ProcessorButton getFatherAnchor() {return (ProcessorButton) super.getFatherAnchor();}

	@Override
	protected ProcessorButton getSonAnchor() {return (ProcessorButton) super.getSonAnchor();}

	
	@Override
	protected void finalize() {
		panel.remove(this);
		getFatherAnchor().removeComponentListener(this);
		getSonAnchor().removeComponentListener(this);
	}

	protected boolean shouldChangeBeRecorded() {return true;}
}
