/***************************************************************************
                          wirelabel.cpp  -  description
                             -------------------
    begin                : Sun February 29 2004
    copyright            : (C) 2004 by Michael Margraf
    email                : michael.margraf@alumni.tu-berlin.de
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "wirelabel.h"
#include "viewpainter.h"
#include "qucs_app.h"
#include "some_font_stuff.h"

#include <QString>
#include <QPainter>

// BUG
#define isHWireLabel       0x4020
#define isVWireLabel       0x4040
#define isNodeLabel        0x4080

WireLabel::WireLabel(const QString& _Name, int, int,
                     int _x1, int _y1, int) : Element()
{
	incomplete();
  // _cx = cx;
  // _cy = cy;
  x1 = _x1;
  y1 = _y1;
  setName(_Name); //?!
  set_label(_Name.toStdString());
  initValue = "";

  // Type = _Type;
  isHighlighted = false;
}

WireLabel::~WireLabel()
{
}

// ----------------------------------------------------------------
void WireLabel::setCenter(int, int, bool)
{
	incomplete();
#if 0
  switch(Type) {
    case isMovingLabel:
      if(relative) {
        x1 += x_;  cx += x_;
        y1 += y_;  cy += y_;
      }
      else {
        x1 = x_;  cx = x_;
        y1 = y_;  cy = y_;
      }
      break;
    case isHMovingLabel:
      if(relative) { x1 += x_;  cx += x_; }
      else { x1 = x_;  cx = x_; }
      break;
    case isVMovingLabel:
      if(relative) { y1 += y_;  cy += y_; }
      else { y1 = y_;  cy = y_; }
      break;
    default:
      if(relative) {
        x1 += x_;  y1 += y_; // moving cx/cy is done by owner (wire, node)
      }
      else { x1 = x_; y1 = y_; }
  }
#endif
}

// ----------------------------------------------------------------
void WireLabel::paint(ViewPainter *p) const
{
  QFont f = p->font(); // save current font
  QFont newFont = f;

  if (isHighlighted) {
//    QColor highlightfill (Qt::blue);
//    highlightfill.setAlpha(50);
//    p->fillRect(x1-1, y1-1, x2, y2, highlightfill);
    p->setPen(QPen(Qt::darkBlue,3));
    newFont.setWeight (QFont::Bold);
  } else {
    newFont.setWeight (QFont::Normal);
    p->setPen(QPen(Qt::black,1));
  }
  p->setFont (newFont);
  int x2=0; int y2=0;
  x2 = p->drawText("TODO", x1, y1, &y2);
  p->setFont(f); // restore old font

  int xpaint=0, ypaint=4, phi=0;
//  switch(Type) {
//    case isVWireLabel:  ypaint=0; xpaint=4; phi=16*140; break;
//    case isHWireLabel:  phi=16*50; break;
//    case isNodeLabel:   ypaint = 0;
//    default:            ;
//  }

  int c, d;
  int a = int(double(x2) / p->Scale) >> 1;
  int b = int(double(y2) / p->Scale) >> 1;
  if(cx() < x1+a) {    // where should frame be painted ?
    if(cy() < y1+b) {
      if(phi == 16*50)  phi += 16*180;
      p->map(x1-3, y1-2, a, b);    // low right
      c = a + (x2>>1);
      d = b + y2;
      p->map(cx()+xpaint, cy()+ypaint, xpaint, ypaint);
    }
    else {
      if(phi != 0)  phi += 16*180;
      p->map(x1-3, y1+1, a, b);    // up right
      b += y2;
      c  = a + (x2>>1);
      d  = b - y2;
      p->map(cx()+xpaint, cy()-ypaint, xpaint, ypaint);
    }
  }
  else {
    if(cy() < y1+b) {
      p->map(x1+3, y1-2, a, b);   // low left
      a += x2;
      c  = a - (x2>>1);
      d  = b + y2;
      p->map(cx()-xpaint, cy()+ypaint, xpaint, ypaint);
    }
    else {
      if(phi > 16*90)  phi += 16*180;
      p->map(x1+3, y1+1, a, b);    // up left
      a += x2;
      b += y2;
      c  = a - (x2>>1);
      d  = b - y2;
      p->map(cx()-xpaint, cy()-ypaint, xpaint, ypaint);
    }
  }

  if(initValue.isEmpty())
    p->setPen(QPen(Qt::darkMagenta,0));
  else
    p->setPen(QPen(Qt::red,0));

  if(phi)  p->drawArc(cx()-4, cy()-4, 8, 8, phi, 16*255);
  p->drawLine(a, b, c, b);
  p->drawLine(a, b, a, d);
  p->drawLine(xpaint, ypaint, a, b);

  x2 = int(double(x2) / p->Scale);
  y2 = int(double(y2) / p->Scale);

#if 0
  if(isSelected()) {
    p->setPen(QPen(Qt::darkGray,3));
    p->drawRoundRect(x1-2, y1-2, x2+6, y2+5);
  }
#endif
}

// ----------------------------------------------------------------
void WireLabel::setName(const QString& Name_)
{
  //setTypeName("wirelabel"); //  = Name_; //?!
  set_label(Name_.toStdString());
  
  // get size of text using the screen-compatible metric
//  FontMetrics metrics;
//  QSize r = metrics.size(0, Name_);
  // x2 = r.width();
  // y2 = r.height()-2;    // remember size of text
}

// ----------------------------------------------------------------
// Converts all necessary data of the wire into a string. This can be used to
// save it to an ASCII file or to transport it via the clipboard.
// Wire labels use the same format like wires, but with length zero.
#if 0
QString WireLabel::save()
{
	unreachable();
  QString s("<");
	s += QString::number(cx())+" "+QString::number(cy())+" "
	  +  QString::number(cx())+" "+QString::number(cy())
	  +  " \""+label() +"\" "
	  +  QString::number(x1)+" "+QString::number(y1)+" 0 \""
	  +  initValue+"\">";
  return s;
}
#endif


void WireLabel::getLabelBounding(int& , int&, int& , int& )
{
	incomplete();
//    _xmin = std::min(x1,x1+(x2+6));
//    _xmax = std::max(x1,x1+(x2+6));
//    _ymin = std::min(y1,y1+(y2+6));
//    _ymax = std::max(y1,y1+(y2+5));
//    _ymax = std::max(cy(),_ymax);
}