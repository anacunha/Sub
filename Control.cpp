// Control.cpp : implementation file
//

#include "stdafx.h"
#include "Sub.h"
#include "Control.h"
#include "SubDoc.h"
#include "SubView.h"

#include "average.hh"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

/////////////////////////////////////////////////////////////////////////////
// CControl dialog


CControl::CControl( CSubView* pParent /*=NULL*/)
	: CDialog(CControl::IDD, (CWnd*)pParent), view( pParent )
{
	//{{AFX_DATA_INIT(CControl)
	m_Plane = TRUE;
	m_Edges = FALSE;
	m_Vertices = TRUE;
	m_EditLoRes = FALSE;
	m_MoveLoRes = FALSE;
	m_SurfaceFilled = TRUE;
	m_Wireframe = FALSE;
	//}}AFX_DATA_INIT
}


void CControl::DoDataExchange(CDataExchange* pDX)
{
	BOOL old_SurfaceFilled = m_SurfaceFilled;
	BOOL old_Wireframe = m_Wireframe;
	BOOL old_Plane = m_Plane;
	BOOL old_Edges = m_Edges;
	BOOL old_Vertices = m_Vertices;
	CDialog::DoDataExchange(pDX);
	//{{AFX_DATA_MAP(CControl)
	DDX_Control(pDX, IDC_SUBSTYLE, m_Substyle);
	DDX_Control(pDX, IDC_LEVELSLIDER, m_LevelSlider);
	DDX_Control(pDX, IDC_LIMIT, m_Limit);
	DDX_Check(pDX, IDC_PLANE, m_Plane);
	DDX_Check(pDX, IDC_EDGES, m_Edges);
	DDX_Check(pDX, IDC_VERTICES, m_Vertices);
	DDX_Check(pDX, IDC_EDITLORES, m_EditLoRes);
	DDX_Check(pDX, IDC_MOVELORES, m_MoveLoRes);
	DDX_Check(pDX, IDC_SURFACE_FILLED, m_SurfaceFilled);
	DDX_Check(pDX, IDC_WIREFRAME, m_Wireframe);
	//}}AFX_DATA_MAP
	if(m_Wireframe != old_Wireframe || m_Plane != old_Plane || m_Vertices != old_Vertices || m_SurfaceFilled != old_SurfaceFilled || m_Edges != old_Edges)
		view->Invalidate(FALSE);
}


BEGIN_MESSAGE_MAP(CControl, CDialog)
	//{{AFX_MSG_MAP(CControl)
	ON_WM_HSCROLL()
	ON_CBN_SELCHANGE(IDC_SUBSTYLE, OnSelchangeSubstyle)
	ON_BN_CLICKED(IDC_LIMIT, OnBnClickedLimit)
	//}}AFX_MSG_MAP	
END_MESSAGE_MAP()

/////////////////////////////////////////////////////////////////////////////
// CControl message handlers

void CControl::OnHScroll(UINT nSBCode, UINT nPos, CScrollBar* pScrollBar) 
{
	if( ((CSliderCtrl*)pScrollBar) == &m_LevelSlider && nSBCode != SB_ENDSCROLL )
	{
		view->GetDocument()->SetLevel( ((CSliderCtrl*)pScrollBar)->GetPos(), m_Limit.GetCheck() );
	}
}

BOOL CControl::OnInitDialog() 
{
	CDialog::OnInitDialog();
	
	m_LevelSlider.SetRange( 0, 5 );
	m_Substyle.SetCurSel( 0 );

	return TRUE;
}

BOOL CControl::OnCmdMsg(UINT nID, int nCode, void* pExtra, AFX_CMDHANDLERINFO* pHandlerInfo) 
{
	BOOL val = CDialog::OnCmdMsg(nID, nCode, pExtra, pHandlerInfo);
	// Grab data from controls every time we have a command msg.  This way we get that
	// lovely non-modal effect.
	UpdateData( TRUE );
	return val;
}

void CControl::OnSelchangeSubstyle() 
{
	static int oldstyle = 0;
	int style = m_Substyle.GetCurSel();
	if( style == oldstyle )
		return;

	oldstyle = style;

	switch( style )
	{
		case 0:
			m_Limit.SetCheck(false);
			m_Limit.EnableWindow(false);
			view->GetDocument()->SetAverage( new AvgNOOP() );
			break;
		case 1:
			m_Limit.SetCheck(false);
			m_Limit.EnableWindow(false);
			view->GetDocument()->SetAverage( new AvgAdHoc( false ) );
			break;
		case 2:
			m_Limit.SetCheck(false);
			m_Limit.EnableWindow(false);
			view->GetDocument()->SetAverage( new AvgAdHoc( true ) );
			break;
		case 3:
			OutputDebugString("Loop\n");
			m_Limit.EnableWindow(true);
			view->GetDocument()->SetAverage( new AvgEval( false ) );
			break;
		case 4:
			OutputDebugString("Butterfly\n");
			m_Limit.SetCheck(false);
			m_Limit.EnableWindow(false);
			view->GetDocument()->SetAverage( new AvgEval( true ) );
			break;
	}
}

void CControl::OnOK()
{
}

void CControl::OnBnClickedLimit()
{
	view->GetDocument()->PushToLimit(m_Limit.GetCheck());
}
