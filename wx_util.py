import wx

def wxButton2int( wx_mouse_button ):
    if wx_mouse_button == wx.MOUSE_BTN_LEFT: return 0
    elif wx_mouse_button == wx.MOUSE_BTN_RIGHT: return 1
    elif wx_mouse_button == wx.MOUSE_BTN_MIDDLE: return 2
    else: return -1

def GetMousePosition( window ):
    return window.ScreenToClient( wx.GetMousePosition() )

def SaveWindowPositionAndSize( frame, window_name ):
    prefs = wx.ConfigBase.Get()
    
    pos = list(map( str, frame.GetPosition() ))
    size = list(map( str, frame.GetSize() ))
    prefs.Write( "/windows/" + window_name + "/x", pos[0] )
    prefs.Write( "/windows/" + window_name + "/y", pos[1] )
    prefs.Write( "/windows/" + window_name + "/w", size[0] )
    prefs.Write( "/windows/" + window_name + "/h", size[1] )
    
    prefs.Flush()

def LoadWindowPositionAndSize( window_name ):
    prefs = wx.ConfigBase.Get()
    
    pos = prefs.Read( "/windows/" + window_name + "/x" ), prefs.Read( "/windows/" + window_name + "/y" )
    size = prefs.Read( "/windows/" + window_name + "/w" ), prefs.Read( "/windows/" + window_name + "/h" )
    
    pos = list(map( int, pos ))
    size = list(map( int, size ))
    
    screen_width = wx.SystemSettings.GetMetric( wx.SYS_SCREEN_X )
    screen_height = wx.SystemSettings.GetMetric( wx.SYS_SCREEN_Y )
    
    size[0] = min( size[0], screen_width )
    size[1] = min( size[1], screen_height )
    pos[0] = max( 0, min( pos[0], screen_width - size[0] ) )
    pos[1] = max( 0, min( pos[1], screen_height - size[1] ) )
    
    return pos, size

def SaveDirectoryPath( last_wd, use_name ):
    prefs = wx.ConfigBase.Get()
    
    prefs.Write( "/directories/" + use_name, last_wd )
    
    prefs.Flush()

def LoadDirectoryPath( use_name ):
    prefs = wx.ConfigBase.Get()
    
    last_wd = prefs.Read( "/directories/" + use_name )
    
    if last_wd == '': raise ValueError
    return last_wd

class _Timers( object ):
    def __init__( self ):
        self.__timers = {}
        self.__timerNextID = 1
    
    def RequestTimer( self, callback, timeout_millis, once = True ):
        ## My timer callback class
        class TimerCB( wx.Timer ):
            def __init__( self, callback ):
                wx.Timer.__init__( self )
                self.callback = callback
            
            def Notify( self ):
                self.callback()
        
        timer_id = self.__timerNextID
        
        timer = TimerCB( callback )
        timer.Start( timeout_millis, once )
        
        self.__timers[ timer_id ] = timer
        self.__timerNextID += 1
        
        return timer_id
    
    def CancelTimer( self, timer_id ):
        self.__timers[ timer_id ].Stop()
        del self.__timers[ timer_id ]

## This is an easy way to implement a singleton in Python.
## It can still be subclassed:
## class Foo( Timers.__class__ ):
Timers = _Timers()
del _Timers
