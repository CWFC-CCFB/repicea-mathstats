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
import java.awt.Container;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.Window;

import repicea.gui.REpiceaShowableUIWithParent;
import repicea.gui.UIControlManager;

@SuppressWarnings("serial")
public class ProcessorLinkLine extends ValidProcessorLinkLine implements REpiceaShowableUIWithParent {

	private final SystemLabel label;
	private transient ProcessorLinkLineSlider guiInterface;
	
	protected ProcessorLinkLine(SystemPanel panel, Processor fatherProcessor, Processor sonProcessor) {
		super(panel, fatherProcessor.getUI(panel), sonProcessor.getUI(panel));
		label = new SystemLabel();
		fatherProcessor.addSubProcessor(sonProcessor);
		fatherProcessor.getUI(panel).addComponentListener(this);
		sonProcessor.getUI(panel).addComponentListener(this);
		addMouseListener(new SystemComponentMouseAdapter(this));
		panel.add(label);
		setLabel();
	}
	
	@Override
	protected void draw(Graphics g) {
		Processor fatherProcessor = getFatherAnchor().getOwner();
		Processor sonProcessor = getSonAnchor().getOwner();
		if (!fatherProcessor.isValid()) {
			g.setColor(Color.RED);
		} else if (fatherProcessor.isPartOfEndlessLoop() && sonProcessor.isPartOfEndlessLoop()) {
			g.setColor(Color.ORANGE);
		} else {
			g.setColor(Color.BLACK);
		}
		super.draw(g);
		g.setColor(UISetup.DefaultColor);
		Point point = getLocation();
		label.setLocation(new Point(point.x + getSize().width, point.y + getSize().height));
	}
	
	@Override
	protected void finalize() {
		super.finalize();
		panel.remove(label);
		getFatherAnchor().getOwner().removeSubProcessor(getSonAnchor().getOwner());
	}
	
	protected void setLabel() {
		Processor fatherProcessor = getFatherAnchor().getOwner();
		Processor sonProcessor = getSonAnchor().getOwner();
		if (fatherProcessor.getSubProcessorIntakes().get(sonProcessor) == null) {
			fatherProcessor.getSubProcessorIntakes().get(sonProcessor);
		}
		int intake = fatherProcessor.getSubProcessorIntakes().get(sonProcessor);
		label.setText(intake + "%");
		panel.repaint();
	}

	@Override
	public ProcessorLinkLineSlider getUI(Container container) {
		if (guiInterface == null) {
			guiInterface = new ProcessorLinkLineSlider((Window) container, this);
		}
		return guiInterface;
	}

	@Override
	public void showUI(Window window) {
		Point referencePoint = getLocationOnScreen();
		ProcessorLinkLineSlider dlg = getUI(window);
		Point currentLocation = new Point(referencePoint.x + 20, referencePoint.y - dlg.getSize().height - 20);
		UIControlManager.setLocation(dlg, currentLocation);
		dlg.setVisible(true);
	}
	
}
