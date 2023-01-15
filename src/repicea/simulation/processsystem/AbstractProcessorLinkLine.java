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

import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Point;
import java.awt.geom.Path2D;

@SuppressWarnings("serial")
public abstract class AbstractProcessorLinkLine extends AbstractPermissionProviderButton {

	private static final int OriginalArrowSize = 15;
	private static final double ArrowRatio = .5;
	private static final double ArrowWaisting = .3;
	
	protected static enum AnchorPosition {
		TOP,
		RIGHT,
		LEFT,
		BOTTOM;
	}
	
	protected final SystemPanel panel;
	private final AnchorProvider fatherAnchor;
	private AnchorProvider sonAnchor;
	
	private static final int BezierFactor = 50;

	
	private AnchorPosition fatherAnchorPosition;
	private AnchorPosition sonAnchorPosition;
	
	protected AbstractProcessorLinkLine(SystemPanel panel, AnchorProvider fatherAnchor, AnchorProvider sonAnchor) {
		super(panel.getListManager().getGUIPermission());
		this.panel = panel;
		this.fatherAnchor = fatherAnchor;
		setSonAnchor(sonAnchor);
		setAnchorPositions(AnchorPosition.RIGHT, AnchorPosition.LEFT);
	}

	protected void setAnchorPositions(AnchorPosition fatherAnchorPosition, AnchorPosition sonAnchorPosition) {
		this.fatherAnchorPosition = fatherAnchorPosition;
		this.sonAnchorPosition = sonAnchorPosition;
	}

	private final Point getAnchor(AnchorProvider provider, AnchorPosition position) {
		switch(position) {
		case RIGHT:
			return provider.getRightAnchor();
		case LEFT:
			return provider.getLeftAnchor();
		case TOP:
			return provider.getTopAnchor();
		case BOTTOM:
			return provider.getBottomAnchor();
		}
		return null;
	}
	
	protected void setSonAnchor(AnchorProvider sonAnchor) {this.sonAnchor = sonAnchor;}
	
	
	protected void setStroke(Graphics2D g2) {
		if (isSelected()) {
			g2.setStroke(UISetup.BoldStroke);
		} else {
			g2.setStroke(UISetup.DefaultStroke);
		}
	}
	
	protected void draw(Graphics g) {
		Graphics2D g2 = (Graphics2D) g;
		Point fatherLocation = getAnchor(fatherAnchor, fatherAnchorPosition);
		Point sonLocation = getAnchor(sonAnchor, sonAnchorPosition);
		if (sonLocation != null) {
			int midX = (int) ((fatherLocation.x + sonLocation.x) * .5);
			int midY = (int) ((fatherLocation.y + sonLocation.y) * .5);
			setLocation(midX - getSize().width / 2, midY - getSize().height / 2);
			setStroke(g2);
			int controlX1 = midX;
			int controlY1 = fatherLocation.y;
			int controlX2 = midX;
			int controlY2 = sonLocation.y;
			if (sonLocation.x - fatherLocation.x < SystemLayout.convertToRelative(BezierFactor)) {
				controlX1 = fatherLocation.x + SystemLayout.convertToRelative(BezierFactor);
				controlX2 = sonLocation.x - SystemLayout.convertToRelative(BezierFactor);
				int offset = SystemLayout.convertToRelative(BezierFactor);
				if (fatherLocation.y > sonLocation.y) {
					offset = -SystemLayout.convertToRelative(BezierFactor);
				}
				controlY1 = fatherLocation.y + offset;
				controlY2 = sonLocation.y - offset;
			}
			Path2D.Double curve = new Path2D.Double();
			curve.moveTo(fatherLocation.x, fatherLocation.y);
			curve.curveTo(controlX1, controlY1,
					controlX2, controlY2,
					sonLocation.x, sonLocation.y);
			g2.draw(curve);

			double x95 = getBezierCoordinate(controlX1, controlX2, sonLocation.x, 0.90);
			double y95 = getBezierCoordinate(controlY1, controlY2, sonLocation.y, 0.90);
			double angle = Math.atan2(sonLocation.y - y95, x95 - sonLocation.x);

//            double arrowRatio = 0.5;
            double arrowLength = SystemLayout.convertToRelative(OriginalArrowSize);

            Path2D.Double arrow = new Path2D.Double();

            double waisting = ArrowWaisting * arrowLength;
            
            arrow.moveTo(sonLocation.x, sonLocation.y);
            lineTo(arrow, sonLocation.x - arrowLength, 
            		sonLocation.y + arrowLength * ArrowRatio, 
            		sonLocation.x, sonLocation.y, angle);
            quadTo(arrow, sonLocation.x - arrowLength + waisting, sonLocation.y,
            		sonLocation.x - arrowLength, sonLocation.y - arrowLength * ArrowRatio,
            		sonLocation.x, sonLocation.y,  angle);
            arrow.lineTo (sonLocation.x, sonLocation.y);
 
            g2.fill (arrow);

			g2.setStroke(UISetup.DefaultStroke);		// to the default value
        }
	}

	protected AnchorProvider getFatherAnchor() {return fatherAnchor;}
	protected AnchorProvider getSonAnchor() {return sonAnchor;}
	
	protected boolean contains(Object button) {
		return getFatherAnchor().equals(button) || getSonAnchor().equals(button);
	}
	
	@Override
	protected void finalize() {}

	@Override
	public boolean equals(Object obj) {
		if (obj instanceof AbstractProcessorLinkLine) {
			AbstractProcessorLinkLine thatLinkLine = (AbstractProcessorLinkLine) obj;
			if (getFatherAnchor().equals(thatLinkLine.getFatherAnchor())) {
				if (getSonAnchor().equals(thatLinkLine.getSonAnchor())) {
					return true;
				}
			}
		}
		return false;
	}

	
	
	double getBezierCoordinate(double x0, double x1, double x2, double t) {
		double x = x0 * (1 - t) * (1 - t) + 2 * x1 * t * (1 - t) + x2 * t * t;
		return x;
	}

	double[] rotatePointAroundReference(double x, double y, double xRef, double yRef, double angle) {
		double xDiff = x - xRef;
		double yDiff = y - yRef;
		return new double[] {xDiff * Math.cos(angle - Math.PI) + yDiff * Math.sin(angle - Math.PI) + xRef, 
				-xDiff * Math.sin(angle - Math.PI) + yDiff * Math.cos(angle - Math.PI) + yRef} ;
	}
	
	void lineTo(Path2D.Double arrow, double x, double y, double xRef, double yRef, double angle) {
		double[] newCoordinate = rotatePointAroundReference(x, y, xRef, yRef, angle);
		arrow.lineTo(newCoordinate[0], newCoordinate[1]);
	}

	void quadTo(Path2D.Double arrow, double x1, double y1, double x2, double y2, double xRef, double yRef, double angle) {
		double[] newCoordinate1 = rotatePointAroundReference(x1, y1, xRef, yRef, angle);
		double[] newCoordinate2 = rotatePointAroundReference(x2, y2, xRef, yRef, angle);
		arrow.quadTo(newCoordinate1[0], newCoordinate1[1], newCoordinate2[0], newCoordinate2[1]);
	}

}
