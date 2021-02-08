/***************************************************************************
                              graphictext.cpp
                             -----------------
    copyright            : (C) 2003 by Michael Margraf
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 3 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include "qucs_app.h"
#include "mnemo.h"
#include "viewpainter.h"
#include "graphictext.h"
#include "graphictextdialog.h"
#include "schematic_doc.h" // BUG
#include "misc.h"
#include "qucs_globals.h"
#include "module.h"
#include "some_font_stuff.h"

#include <QPainter>
#include <QPushButton>
#include <QLineEdit>
#include <QTextEdit>

#include <cmath>

#include "../legacy_painting.h"

#include <QColor>
#include <QFont>
#include <QString>


namespace{
class GraphicText : public LegacyPainting {
public:
  GraphicText();
 ~GraphicText();

  void paintScheme(SchematicDoc*);
  void getCenter(int&, int&);
//  void setCenter(int, int, bool relative=false);

  Element* newOne() const {return new GraphicText(*this);}
  Element* clone() const {return new GraphicText(*this);}

  static Element* info(QString&, char* &, bool getNewOne=false);
  bool load(const QString&);
  QString save();
  QString saveCpp();
  QString saveJSON();
  void paint(ViewPainter*);
  void MouseMoving(SchematicDoc*, int, int, int, int, SchematicDoc*, int, int, bool);
  bool MousePressing();
  bool getSelected(float, float, float);
  void Bounding(int&, int&, int&, int&);

  void rotate();
  void mirrorX();
  void mirrorY();
  bool Dialog();

  QColor   Color;
  QFont    Font;
  QString  Text;
  int      Angle;
};

GraphicText D;
Dispatcher<Element>::INSTALL p(&element_dispatcher, "Text", &D);
Module::INSTALL pp("paintings", &D);
}

GraphicText::GraphicText() : LegacyPainting()
{
  Name = "Text ";
  Color = QColor(0,0,0);
  Font = QString_(QucsSettings.font);
  x1 = x2 = 0;
  y1 = y2 = 0;
  Angle = 0;
}

GraphicText::~GraphicText()
{
}

// -----------------------------------------------------------------------
void GraphicText::paint(ViewPainter *p)
{

	// return; ?
//	   auto cx=Element::cx();
//     auto cy=Element::cy();

  // keep track of painter state
  p->save();

  auto wm = p->worldMatrix();
#if 0
  QMatrix Mat(1.0, 0.0, 0.0, 1.0, p->DX + float(cx) * p->Scale,
				   p->DY + float(cy) * p->Scale);
  p->setWorldMatrix(Mat);
  p->rotate(-Angle);   // automatically enables transformation

#endif
  int Size = Font.pointSize();
  Font.setPointSizeF( p->FontScale * float(Size) );

  QFont f = p->font();
#if 0
  p->setPen(Color);
  p->setFont(Font);

  // Because of a bug in Qt 3.1, drawing this text is dangerous, if it
  // contains linefeeds. Qt has problems with linefeeds. It remembers the
  // last font metrics (within the font ???) and does not calculate it again.
  // The error often appears at a very different drawText function !!!
#endif
  int w, h;
  w = p->drawTextMapped(Text, 0, 0, &h);

#if 0
  if(isSelected()) {
    p->setPen(QPen(Qt::darkGray,3));
    p->drawRect(-3, -2, w+6, h+5);
  }
#endif

  Font.setPointSize(Size);   // restore real font size
  p->setWorldMatrix(wm);
  p->setWorldMatrixEnabled(false);

  // restore painter state
  p->restore();

  x2 = int(float(w) / p->Scale);
  y2 = int(float(h) / p->Scale);
  p->setFont(f);
}

// -----------------------------------------------------------------------
void GraphicText::paintScheme(SchematicDoc *p)
{
	incomplete();
  // FIXME #warning QMatrix wm = p->worldMatrix();
  // FIXME #warning QMatrix Mat (wm.m11(), 0.0, 0.0, wm.m22(),
// FIXME #warning 		wm.dx() + double(cx) * wm.m11(),
// FIXME #warning 		wm.dy() + double(cy) * wm.m22());
  // FIXME #warning p->setWorldMatrix(Mat);
  // FIXME #warning p->rotate(-Angle);
  p->PostPaintEvent(_Rect, 0, 0, x2, y2);

  // FIXME #warning p->setWorldMatrix(wm);
}

// ------------------------------------------------------------------------
void GraphicText::getCenter(int& x, int &y)
{
	   auto cx=Element::cx();
     auto cy=Element::cy();

  x = cx+(x2>>1);
  y = cy+(y2>>1);
}

// -----------------------------------------------------------------------
// Sets the center of the painting to x/y.
// void GraphicText::setCenter(int x, int y, bool relative)
// {
//   if(relative) {
// 	  _cx += x;
// 	  _cy += y;
//   }else{
// 	  _cx = x-(x2>>1);
// 	  _cy = y-(y2>>1);
//   }
// }
// 
// -----------------------------------------------------------------------
Element* GraphicText::info(QString& Name, char* &BitmapFile, bool getNewOne)
{
  Name = QObject::tr("Text");
  BitmapFile = (char *) "text";

  if(getNewOne)  return new GraphicText();
  return 0;
}

// -----------------------------------------------------------------------
bool GraphicText::load(const QString& s)
{
  bool ok;

  QString n;
  n  = s.section(' ',1,1);    // cx
  int cx = n.toInt(&ok);
  if(!ok) return false;

  n  = s.section(' ',2,2);    // cy
  int cy = n.toInt(&ok);
  if(!ok) return false;

  setCenter(pos_t(cx, cy));

  n  = s.section(' ',3,3);    // Size
  Font.setPointSize(n.toInt(&ok));
  if(!ok) return false;

  n  = s.section(' ',4,4);    // Color
  Color.setNamedColor(n);
  if(!Color.isValid()) return false;

  n  = s.section(' ',5,5);    // Angle
  Angle = n.toInt(&ok);
  if(!ok) return false;

  Text = s.mid(s.indexOf('"')+1);    // Text (can contain " !!!)
  Text.truncate(Text.length()-1);
  if(Text.isEmpty()) return false;

  misc::convert2Unicode(Text);
  // get size of text using the screen-compatible metric
  FontMetrics metrics;
  QSize r = metrics.size(0, Text);    // get overall size of text
  x2 = r.width();
  y2 = r.height();

  return true;
}

// -----------------------------------------------------------------------
QString GraphicText::save()
{
	   auto cx=Element::cx();
     auto cy=Element::cy();

  QString t = Text;
  misc::convert2ASCII(t);

  // The 'Text' property has to be the last within the line !
  QString s = Name+QString::number(cx)+" "+QString::number(cy)+" "
		+ QString::number(Font.pointSize())+" "+Color.name()+" "
		+ QString::number(Angle) + " \""+t+"\"";
  return s;
}

// --------------------------------------------------------------------------
QString GraphicText::saveCpp()
{
	   auto cx=Element::cx();
     auto cy=Element::cy();

  QString t = Text;
  misc::convert2ASCII(t);

  QString s =
    QString ("new Text (%1, %2, \"%3\", QColor (\"%4\"), %5, %6, %7)").
    arg(cx).arg(cy).arg(t).
    arg(Color.name()).arg(Font.pointSize()).
    arg(cos(pi * Angle / 180.0)).arg(sin(pi * Angle / 180.0));
  s = "Texts.append (" + s + ");";
  return s;
}

QString GraphicText::saveJSON()
{
	   auto cx=Element::cx();
     auto cy=Element::cy();

  QString t = Text;
  misc::convert2ASCII(t);

  QString s =
    QString ("{\"type\" : \"graphictext\", "
      "\"x\" : %1, \"y\" : %2, \"s\" : \"%3\", "
      "\"color\" : \"%4\", \"size\" : %5, \"cos\" : %6, \"sin\" : %7},").
      arg(cx).arg(cy).arg(t).
      arg(Color.name()).arg(Font.pointSize()).
      arg(cos(pi * Angle / 180.0)).arg(sin(pi * Angle / 180.0));
  return s;
}

// -----------------------------------------------------------------------
// fx/fy are the precise coordinates, gx/gy are the coordinates set on grid.
// x/y are coordinates without scaling.
void GraphicText::MouseMoving(
	SchematicDoc*, int, int, int , int,
	SchematicDoc *, int , int, bool)
{
#if 0
  // FIXME #warning p->setPen(Qt::SolidLine);
  if(drawn) {
    p->PostPaintEvent(_Line, x1+15, y1+15, x1+20, y1,0,0,true);  // erase old cursor symbol
    p->PostPaintEvent(_Line, x1+26, y1+15, x1+21, y1,0,0,true);
    p->PostPaintEvent(_Line, x1+17, y1+8,  x1+23, y1+8,0,0,true);
  }
  x1 = x;
  y1 = y;
  p->PostPaintEvent(_Line, x1+15, y1+15, x1+20, y1,0,0,true);  // paint new cursor symbol
  p->PostPaintEvent(_Line, x1+26, y1+15, x1+21, y1,0,0,true);
  p->PostPaintEvent(_Line, x1+17, y1+8,  x1+23, y1+8,0,0,true);

  cx = gx;
  cy = gy;
#endif
}

// ------------------------------------------------------------------------
bool GraphicText::MousePressing()
{
  return Dialog();
}

// ------------------------------------------------------------------------
// Checks if the coordinates x/y point to the painting.
// 5 is the precision the user must point onto the painting.
bool GraphicText::getSelected(float fX, float fY, float)
{
	   auto cx=Element::cx();
     auto cy=Element::cy();

  double phi  = pi/180.0*double(Angle);
  float  sine = sin(phi), cosine = cos(phi);

  fX -= float(cx);
  fY -= float(cy);
  int _x = int( fX*cosine - fY*sine );
  int _y = int( fY*cosine + fX*sine );

  if(_x >= 0) if(_y >= 0) if(_x <= x2) if(_y <= y2)
    return true;

  return false;
}

// ------------------------------------------------------------------------
void GraphicText::Bounding(int& xmin, int& ymin, int& xmax, int& ymax)
{
	   auto cx=0;
     auto cy=0;

  double phi = pi/180.0*double(Angle);
  double sine = sin(phi), cosine = cos(phi);
  int dx = int( double(y2) * sine );
  int dy = int( double(y2) * cosine );
  xmin = dx;  xmax = cx;
  ymin = dy;  ymax = cy;
  if(xmin < 0)  xmin += cx;
  else { xmax += xmin;  xmin = cx; }
  if(ymin < 0)  ymin += cy;
  else { ymax += ymin;  ymin = cy; }

  int x = cx + int( double(x2) * cosine );
  if(xmax < x)  xmax = x;
  else if(xmin > x)  xmin = x;
  x += dx;
  if(xmax < x)  xmax = x;
  else if(xmin > x)  xmin = x;

  int y = cy - int( double(x2) * sine );
  if(ymax < y)  ymax = y;
  else if(ymin > y)  ymin = y;
  y += dy;
  if(ymax < y)  ymax = y;
  else if(ymin > y)  ymin = y;
}

// -----------------------------------------------------------------------
// Rotates around the center.
void GraphicText::rotate()
{
  Angle += 90;
  Angle %= 360;
  // _cx -= x2 >> 1;
  // _cy -= y2 >> 1;
}

// -----------------------------------------------------------------------
// Mirrors about center line.
void GraphicText::mirrorX()
{   // do not mirror, because unreadable
}

// -----------------------------------------------------------------------
// Mirrors about center line.
void GraphicText::mirrorY()
{   // do not mirror, because unreadable
}

// -----------------------------------------------------------------------
// Calls the property dialog for the painting and changes them accordingly.
// If there were changes, it returns 'true'.
//
// schematicWidget?
bool GraphicText::Dialog()
{
  QFont f(QString_(QucsSettings.font));   // to avoid wrong text width
  bool changed = false;

  GraphicTextDialog *d = new GraphicTextDialog();

  // initialize dialog picker color
  misc::setPickerColor(d->ColorButt, Color);

  d->TextSize->setText(QString::number(Font.pointSize()));
  d->Angle->setText(QString::number(Angle));
  QString _Text = Text;
  decode_String(_Text);  // replace special characters with LaTeX commands
  d->text->setText(_Text);

  if(d->exec() == QDialog::Rejected) {
    delete d;
    return false;
  }

  QColor curColor = misc::getWidgetBackgroundColor(d->ColorButt);
  if(Color != curColor) {
    Color = curColor;
    changed = true;
  }
  f.setPointSize(d->TextSize->text().toInt());   // to avoid wrong text width
  if(Font.pointSize()  != d->TextSize->text().toInt()) {
    Font.setPointSize(d->TextSize->text().toInt());
    changed = true;
  }
  int tmp = d->Angle->text().toInt();
  if(Angle  != tmp) {
    Angle = tmp % 360;
    changed = true;
  }

  encode_String(d->text->toPlainText(), _Text);  // create special characters
  if(!_Text.isEmpty())
    if(_Text != Text) {
      Text = _Text;
      changed = true;
    }

  // get font metric using the screen-compatible metric
  QFontMetrics  m(f, 0);
  QSize s = m.size(0, Text); // get size of text
  x2 = s.width();
  y2 = s.height();

  delete d;
  return changed;
}
