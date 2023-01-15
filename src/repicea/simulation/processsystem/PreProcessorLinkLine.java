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
import java.awt.Component;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;

@SuppressWarnings("serial")
public class PreProcessorLinkLine extends AbstractProcessorLinkLine {

	
	protected class InternalMouseMotionListener extends MouseAdapter implements AnchorProvider {
		
		private Point lastKnownLocation;
		private boolean onGoingAction;
		
		@Override
		public Point getRightAnchor() {return lastKnownLocation;}

		@Override
		public Point getLeftAnchor() {return lastKnownLocation;}

		@Override
		public Point getBottomAnchor() {return lastKnownLocation;}

		@Override
		public Point getTopAnchor() {return lastKnownLocation;}

		@Override
		public void mouseDragged(MouseEvent e) {
			onGoingAction = true;
			Component originComponent = (Component) e.getSource();
			if (PreProcessorLinkLine.this.getSonAnchor() instanceof InternalMouseMotionListener) {
				lastKnownLocation = new Point(originComponent.getLocation().x + e.getPoint().x, originComponent.getLocation().y + e.getPoint().y);
			}
			panel.repaint();
		}
		
		@Override
		public void mouseEntered(MouseEvent e) {
			if (onGoingAction) {
				if (e.getSource() instanceof ProcessorButton) {
					ProcessorButton button = (ProcessorButton) e.getSource();
					if (!button.equals(PreProcessorLinkLine.this.getFatherAnchor()) && button.getGUIPermission().isDropGranted()) {
						PreProcessorLinkLine.this.setSonAnchor(button);
					}
				}
			}
		}
		
		@Override
		public void mouseExited(MouseEvent e) {
			if (onGoingAction) {
				PreProcessorLinkLine.this.setSonAnchor(this);
			}
		}

		@Override
		public void mouseReleased(MouseEvent e) {
			onGoingAction = false;
			if (PreProcessorLinkLine.this.getSonAnchor() instanceof ProcessorButton) {
				panel.doLinkAction();
			} else {
				panel.registerLinkBeingCreated(null);
			}
			panel.repaint();
		}
	}
	
	
	private final InternalMouseMotionListener mouseAdapter;
	
	protected PreProcessorLinkLine(SystemPanel panel, AnchorProvider fatherAnchor) {
		super(panel, fatherAnchor, null);
		mouseAdapter = new InternalMouseMotionListener();
		this.setSonAnchor(mouseAdapter);
		for (ProcessorButton processorButton : panel.processorButtons) {
			processorButton.addMouseListener(mouseAdapter);
			processorButton.addMouseMotionListener(mouseAdapter);
		}
	}

	@Override
	protected void finalize() {
		for (ProcessorButton processorButton : panel.processorButtons) {
			processorButton.removeMouseListener(mouseAdapter);
			processorButton.removeMouseMotionListener(mouseAdapter);
		}
	}

	@Override
	protected void setStroke(Graphics2D g2) {
		g2.setColor(Color.BLACK);
		g2.setStroke(UISetup.BoldStroke);
	}

	protected ValidProcessorLinkLine convertIntoProcessorLinkLine() {
		return new ProcessorLinkLine(panel, getFatherAnchor().getOwner(), ((ProcessorButton) getSonAnchor()).getOwner());
	}


	@Override
	protected ProcessorButton getFatherAnchor() {return (ProcessorButton) super.getFatherAnchor();}

}
