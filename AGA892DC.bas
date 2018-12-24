Attribute VB_Name = "AGA892DC"
'函数原型：
'double xipt2Z_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t)
'double xipt2rhom_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t)
'double xipt2rhom_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t)
'
'/*
'AGA892DC:
'    //xi:摩尔分数,p:绝对压力,MPa,t:摄氏度，C
'    /*
'    p:0<p<=12MPa
'    t+273.15:263-338
'    x(CH4):0.7-1.0
'    x(N2):0-0.7
'    x(CO2):0-0.2
'    x(C2H6):0-0.1
'    x(C3H8):0-0.035
'    x(C4H10):0-0.015
'    x(C5H12):0-0.005
'    x(C6H14):0-0.001
'    x(C7H16):0-0.0005
'    x(C8+):0-0.0005
'    x(H2):0-0.1
'    x(CO):0-0.03
'    x(He):0-0.005
'    x(H2O):0-0.00015
'    x(O2):0-0.00015?
'    x(H2S):0-0.00015?
'    x(Ar+Ne+Kr+Xe):0-0.00015?
'    sum(xi)==1.0
'    ID:         1        2  3   4       5       6   7   8   9   10
'    Component:  CH4     N2  CO2 C2H6    C3H8    H2O H2S H2  CO  O2
'    ID:             11      12      13      14      15      16      17  1   8       19      20  21
'   Component:    i-C4H10   n-C4H10 i-C5H12 n-C5H12 C6H14   C7H16   C8H18   C9H20   C10H22  He  Ar
    
'*/
'NOTE: C/C++ int type <==> vba long type ,vba pass array parameter such as:
'Dim ID(1:21)  as Integer
'func(ID(1)) will pass ID array for func
Option Base 1

Public Declare PtrSafe Function xipt2Z_AGA Lib "NGCompressionFactor.dll" Alias "xipt2Z_AGA@28" _
(ByRef xi_raw As Double, ByRef ID_raw As Long, ByVal num_x As Long, ByVal p As Double, ByVal t As Double) As Double

Private Function InputProcess(xi As range, Component As range, xi_raw() As Double, ID_raw() As Long, numx As Long) As Boolean
    'xi Component必须为1维的列/行数组
    Dim Component_raw() As String
    Dim xi0
    Dim Component0
    xi0 = xi
    Component0 = Component
    Dim xi0_flag As Boolean, Component0_flag As Boolean
    xi0_flag = IsArray(xi0)
    Component0_flag = IsArray(Component0)
    If ((xi0_flag = False And Component0_flag = True) Or (xi0_flag = True And Component0_flag = False)) Then
        InputProcess = False
        Exit Function
    ElseIf (xi0_flag = False And Component0_flag = False) Then
        numx = 1
        'ReDim xi_raw(numx)
        ReDim Component_raw(numx)
        xi_raw(1) = xi0
        Component_raw(1) = Component0
    Else
        numx = UBound(xi0, 1)
        If (numx <> UBound(Component0, 1)) Then
            InputProcess = False
            Exit Function
        End If
        If (numx = 1) Then
            numx = UBound(xi0, 2)
            If (numx <> UBound(Component0, 2) Or numx > 21) Then
                InputProcess = False
                Exit Function
            End If
            'ReDim xi_raw(numx)
            ReDim Component_raw(numx)
            For i = 1 To numx
                If (IsNumeric(xi_raw(i)) = False) Then
                    InputProcess = False
                    Exit Function
                End If
                xi_raw(i) = xi0(1, i)
                Component_raw(i) = Component0(1, i)
            Next
        Else
            If (numx > 21) Then
                InputProcess = False
                Exit Function
            End If
            'ReDim xi_raw(numx)
            ReDim Component_raw(numx)
            For i = 1 To numx
                If (IsNumeric(xi_raw(i)) = False) Then
                    InputProcess = False
                    Exit Function
                End If
                xi_raw(i) = xi0(i, 1)
                Component_raw(i) = Component0(i, 1)
            Next
        End If
    End If
    Dim component_list()
    component_list = Array("x(CH4)", "x(N2)", "x(CO2)", "x(C2H6)", "x(C3H8)", "x(H2O)", "x(H2S)", "x(H2)", "x(CO)", "x(O2)", "x(i-C4H10)" _
                         , "x(n-C4H10)", "x(i-C5H12)", "x(n-C5H12)", "x(C6H14)", "x(C7H16)", "x(C8H18)", "x(C9H20)", "x(C10H22)", "x(He)", "x(Ar)")
    On Error Resume Next
    Err.Clear
    'ReDim ID_raw(numx)
    For i = 1 To numx
        ID_raw(i) = Application.WorksheetFunction.IfError(Application.WorksheetFunction.Match(Component_raw(i), component_list, 0), -1)
        If (Err.Number <> 0) Then
            InputProcess = False
            Exit Function
        End If
    Next
    InputProcess = True
End Function

Public Function Z_XiPT_AGA(xi As range, Component As range, p As Double, t As Double)
    'xi 输入的组分摩尔分数的单元格范围，必须为1维列或行
    ' Component 输入的组分名称的单元格范围，形如x(CH4)，必须是与xi对应的列或行
    Dim xi_raw(21) As Double
    Dim ID_raw(21) As Long
    Dim numx As Long
    xi_raw(1) = 0
    ID_raw(1) = 0
    incheck = InputProcess(xi, Component, xi_raw, ID_raw, numx)
    If (incheck = False) Then
        Z_XiPT_AGA = "errorinput"
        Exit Function
    End If
    
    
    Dim strCurDriver As String
    Dim strCurDir As String
    Dim tbkDir As String
    Dim tbkDriver As String
    '备份当前路径
    strCurDir = CurDir
    strCurDriver = Left(strCurDir, 3)
    tbkDir = ThisWorkbook.Path
    tbkDriver = Left(tbkDir, 3)
    '设置工作薄路径为新路径
    ChDrive tbkDriver
    ChDir tbkDir
    
    Z = 0
    Z = xipt2Z_AGA(xi_raw(1), ID_raw(1), numx, p, t)
    
    '恢复之前路径
    ChDrive strCurDriver
    ChDir strCurDir
    Z_XiPT_AGA = Z
End Function






