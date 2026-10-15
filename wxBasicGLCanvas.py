## Based off GLCanvas.py sample code and code from http://wiki.wxpython.org/GLCanvas

import wx
from wx import glcanvas

# The Python OpenGL package can be found at
# http://PyOpenGL.sourceforge.net/
from OpenGL.GL import *
## GLUT launches X11
# from OpenGL.GLUT import glutInit

## Maximum zoom is 16x
kMaxZoom = 16

class BasicGLCanvas( glcanvas.GLCanvas ):
    def __init__(self, parent, pos = wx.DefaultPosition, size = wx.DefaultSize, debugging = True, GLUT = False, **kwargs ):
        attribList = [
            ## RGBA
            glcanvas.WX_GL_RGBA,
            ## Double Buffered
            glcanvas.WX_GL_DOUBLEBUFFER,
            ## 32 bit depth buffer
            ## UPDATE: 32-bits stopped working when I installed the latest wxWidgets with homebrew.
            ##         I have to use 24-bits.
            glcanvas.WX_GL_DEPTH_SIZE, 24,
            ## The list must end with a 0.
            0
            ]
        assert glcanvas.GLCanvas.IsDisplaySupported( attribList )
        
        ## Modern attribute setup: https://wxpython.org/Phoenix/docs/html/wx.glcanvas.GLAttributes.html#wx-glcanvas-glattributes
        ## Only in developer-channel wxPython.
        # dispAttrs = wx.glcanvas.GLAttributes()
        # dispAttrs.PlatformDefaults().MinRGBA(8,8,8,8).DoubleBuffer().Depth(32).EndList()
        
        shared = kwargs.setdefault( 'shared', None )
        #print 'Creating a BasicGLCanvas with the shared context:', shared
        if shared is None:
            glcanvas.GLCanvas.__init__(self, parent, -1, pos = pos, size = size, attribList = attribList )
        else:
            ## This trick comes from the wxPython-users thread "Re: inheritance and GLCanvas with shared context"
            ## accessed at <http://wxpython-users.1045709.n5.nabble.com/inheritance-and-GLCanvas-with-shared-context-td2312431.html>
            self.PostCreate( glcanvas.GLCanvasWithContext( parent, shared, -1, pos = pos, size = size, attribList = attribList ) )
        self.__init = False
        ## self.GetContext() always seems to return None.
        self.context = glcanvas.GLContext( self )
        self.size = None
        self.world_bounds_contains = [ -1, 1, -1, 1 ]
        self.world_bounds = list( self.world_bounds_contains )
        self.debugging = debugging
        
        ## Drawing events
        self.Bind(wx.EVT_ERASE_BACKGROUND, self.__OnEraseBackground)
        self.Bind(wx.EVT_SIZE, self.__OnSize)
        self.Bind(wx.EVT_PAINT, self.__OnPaint)
        
        ## Mouse events
        self.Bind(wx.EVT_LEFT_DOWN, self.__OnMouseDown)
        self.Bind(wx.EVT_RIGHT_DOWN, self.__OnMouseDown)
        self.Bind(wx.EVT_MIDDLE_DOWN, self.__OnMouseDown)
        self.Bind(wx.EVT_LEFT_UP, self.__OnMouseUp)
        self.Bind(wx.EVT_RIGHT_UP, self.__OnMouseUp)
        self.Bind(wx.EVT_MIDDLE_UP, self.__OnMouseUp)
        self.Bind(wx.EVT_MOTION, self.__OnMouseMotion)
        
        ## Keyboard events
        self.Bind(wx.EVT_KEY_DOWN, self.__OnKeyDown)
        self.Bind(wx.EVT_KEY_UP, self.__OnKeyUp)
        
        ## Context menu event
        self.Bind( wx.EVT_CONTEXT_MENU, self.__OnContextMenu )
        
        ## If someone does want GLUT, we need to init it.
        self.GLUT = GLUT
    
    ## private methods:
    
    def __OnEraseBackground(self, event):
        pass # Do nothing, to avoid flashing on Windows.


    def __OnSize(self, event):
        event.Skip()
        
        size = self.GetClientSize()
        #print 'new size:', size
        
        ## Ignore the resize if one of the dimensions is zero.
        if 0 in size: return
        
        self.size = size
        self.__updateWorldBounds()
        
        if self.IsShownOnScreen():
            self.SetCurrent( self.context )
            glViewport(0, 0, size.width, size.height)
    
    def __updateWorldBounds( self ):
        ## These world bounds project squares to squares, and are as small as possible
        ## while containing self.world_bounds_contains.
        
        ## Find the (centered) rect to contain.
        contain = list( self.world_bounds_contains )
        center = ( .5*( contain[0] + contain[1] ), .5*( contain[2] + contain[3] ) )
        contain[0] -= center[0]
        contain[1] -= center[0]
        contain[2] -= center[1]
        contain[3] -= center[1]
        
        contain_width = contain[1] - contain[0]
        contain_height = contain[3] - contain[2]
        
        ## Like parameters to glOrtho: xmin, xmax, ymin, ymax
        self.world_bounds = [ -1, 1, -1, 1 ]
        if self.size.width >= self.size.height:
            aspect = self.size.width / float( self.size.height )
            self.world_bounds[ 0 ] = -aspect
            self.world_bounds[ 1 ] = aspect
            
            contain_aspect = contain_width / contain_height
            if contain_aspect > aspect:
                scale = contain_width / ( 2 * aspect )
            else:
                scale = contain_height / 2
        else:
            aspect = self.size.height / float( self.size.width )
            self.world_bounds[ 2 ] = -aspect
            self.world_bounds[ 3 ] = aspect
            
            contain_aspect = contain_height / contain_width
            if contain_aspect > aspect:
                scale = contain_height / ( 2 * aspect )
            else:
                scale = contain_width / 2
        
        self.world_bounds = [ a * scale for a in self.world_bounds ]
        self.world_bounds[0] += center[0]
        self.world_bounds[1] += center[0]
        self.world_bounds[2] += center[1]
        self.world_bounds[3] += center[1]
        
        self.world_bounds = tuple( self.world_bounds )
        #print "world_bounds_contains:", self.world_bounds_contains
        #print "world_bounds:", self.world_bounds
    
    def GetZoom( self ):
        '''
        untested
        '''
        contain = self.world_bounds_contains
        return 2. / min( contain[1] - contain[0], contain[3] - contain[2] )
    
    def ResetZoom( self ):
        self.world_bounds_contains = [ -1, 1, -1, 1 ]
        
        self.__updateWorldBounds()
    
    def ZoomIn( self, center = None ):
        ## 1 Translate world_bounds_contains to the origin.
        ## 2 Zoom in.
        ## 3 Translate world_bounds_contains to have center 'center'.
        
        ## 1
        contain = list( self.world_bounds_contains )
        cur_center = ( .5*( contain[0] + contain[1] ), .5*( contain[2] + contain[3] ) )
        contain[0] -= cur_center[0]
        contain[1] -= cur_center[0]
        contain[2] -= cur_center[1]
        contain[3] -= cur_center[1]
        
        ## 2
        contain[0] += .25
        contain[1] -= .25
        contain[2] += .25
        contain[3] -= .25
        
        contain[0] = min( -1./kMaxZoom, contain[0] )
        contain[1] = max( 1./kMaxZoom, contain[1] )
        contain[2] = min( -1./kMaxZoom, contain[2] )
        contain[3] = max( 1./kMaxZoom, contain[3] )
        
        ## 3
        if center is None:
            center = cur_center
        
        contain[0] += center[0]
        contain[1] += center[0]
        contain[2] += center[1]
        contain[3] += center[1]
        
        self.world_bounds_contains = contain
        
        self.__updateWorldBounds()
    
    def ZoomOut( self, center = None ):
        ## 1 Translate world_bounds_contains to the origin.
        ## 2 Zoom out.
        ## 3 Translate world_bounds_contains to have center 'center'.
        
        ## 1
        contain = list( self.world_bounds_contains )
        cur_center = ( .5*( contain[0] + contain[1] ), .5*( contain[2] + contain[3] ) )
        contain[0] -= cur_center[0]
        contain[1] -= cur_center[0]
        contain[2] -= cur_center[1]
        contain[3] -= cur_center[1]
        
        ## 2
        contain[0] -= .25
        contain[1] += .25
        contain[2] -= .25
        contain[3] += .25
        
        contain[0] = min( -1./kMaxZoom, contain[0] )
        contain[1] = max( 1./kMaxZoom, contain[1] )
        contain[2] = min( -1./kMaxZoom, contain[2] )
        contain[3] = max( 1./kMaxZoom, contain[3] )
        
        ## 3
        if center is None:
            center = cur_center
        
        contain[0] += center[0]
        contain[1] += center[0]
        contain[2] += center[1]
        contain[3] += center[1]
        
        self.world_bounds_contains = contain
        
        self.__updateWorldBounds()
    
    def ZoomToRect( self, corner_one, corner_two ):
        contain = [
            min( corner_one[0], corner_two[0] ),
            max( corner_one[0], corner_two[0] ),
            min( corner_one[1], corner_two[1] ),
            max( corner_one[1], corner_two[1] )
            ]
        
        ## We'd be done if we didn't want to clamp to the maximum zoom.
        ## Translate to the origin, clamp the zoom, translate back.
        center = ( .5*( contain[0] + contain[1] ), .5*( contain[2] + contain[3] ) )
        contain[0] -= center[0]
        contain[1] -= center[0]
        contain[2] -= center[1]
        contain[3] -= center[1]
        
        contain[0] = min( -1./kMaxZoom, contain[0] )
        contain[1] = max( 1./kMaxZoom, contain[1] )
        contain[2] = min( -1./kMaxZoom, contain[2] )
        contain[3] = max( 1./kMaxZoom, contain[3] )
        
        contain[0] += center[0]
        contain[1] += center[0]
        contain[2] += center[1]
        contain[3] += center[1]
        
        self.world_bounds_contains = contain
        
        self.__updateWorldBounds()
    
    def Pan( self, delta2d ):
        self.world_bounds_contains[0] += delta2d[0]
        self.world_bounds_contains[1] += delta2d[0]
        self.world_bounds_contains[2] += delta2d[1]
        self.world_bounds_contains[3] += delta2d[1]
        
        self.__updateWorldBounds()
    
    def __OnPaint(self, event):
        dc = wx.PaintDC(self)
        self.SetCurrent( self.context )
        if not self.__init:
            self.doInitGL()
            self.init = True
        self.doDraw()
    
    def __OnMouseDown( self, evt ):
        self.CaptureMouse()
        self.doMouseDown( evt )
        ## If the caller didn't want the down event, don't capture the mouse.
        if evt.Skipped: self.ReleaseMouse()
    
    def __OnMouseUp( self, evt ):
        ## Clicks that gave focus to the window will give us mouse up events with no down events.
        if not self.HasCapture():
            evt.Skip()
            return
        self.ReleaseMouse()
        self.doMouseUp( evt )
    
    def __OnMouseMotion(self, evt):
        #if evt.Dragging() and evt.LeftIsDown():
        #    self.lastx, self.lasty = self.x, self.y
        #    self.x, self.y = evt.GetPosition()
        #    self.Refresh(False)
        if evt.Dragging():
            if not self.HasCapture():
                evt.Skip()
                return
            self.doMouseDrag( evt )
        else:
            self.doMouseMove( evt )
    
    def __OnKeyDown( self, evt ):
        self.doKeyDown( evt )
    
    def __OnKeyUp( self, evt ):
        self.doKeyUp( evt )
    
    def __OnContextMenu( self, evt ):
        self.doContextMenu( evt )
    
    ## public methods:
    def GetMousePosition( self ):
        ## NOTE: wx_util.py has a version of this method that takes the
        ##       wxWindow as a parameter.
        return self.ScreenToClient( wx.GetMousePosition() )
    
    def GetMouseState( self ):
        '''
        Returns a wx MouseState object.
        NOTE: The position in the MouseState object will be in screen coordinates.
              It must be passed through self.ScreenToClient() to get event-space coordinates.
        '''
        return wx.GetMouseState()
    
    ## This method is inherited; it requests a redraw of the window.
    #def Refresh()
    
    ## Methods beginning with 'do' are for subclasses to override to
    ## implement functionality.
    def doMouseDown( self, evt ):
        if self.debugging: print('mdown:', evt.GetPosition())
        evt.Skip()
    def doMouseUp( self, evt ):
        if self.debugging: print('mup:', evt.GetPosition())
        evt.Skip()
    def doMouseMove( self, evt ):
        if self.debugging: print('mmove:', evt.GetPosition())
        evt.Skip()
    def doMouseDrag( self, evt ):
        if self.debugging: print('mdrag:', evt.GetPosition())
        evt.Skip()
    
    def doKeyDown( self, evt ):
        if self.debugging: print('kdown:', evt.GetPosition())
        evt.Skip()
    def doKeyUp( self, evt ):
        if self.debugging: print('kup:', evt.GetPosition())
        evt.Skip()
    
    def doContextMenu( self, evt ):
        if self.debugging: print('context menu:', evt.GetPosition())
        evt.Skip()
    
    def doInitGL( self ):
        ## GLUT launches X11
        if self.GLUT:
            import OpenGL.GLUT
            OpenGL.GLUT.glutInit()
        
        glDisable( GL_LIGHTING )
        glDisable( GL_DEPTH_TEST )
        
        glEnable( GL_BLEND )
        glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA )
        
        #glEnable( GL_RESCALE_NORMAL )
        glEnable( GL_NORMALIZE )
        
        #glEnable( GL_POINT_SMOOTH )
        #glEnable( GL_LINE_SMOOTH )
    
    def doDraw( self ):
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(-1, 1, -1, 1, -1, 1)
        
        glMatrixMode(GL_MODELVIEW)
        glLoadIdentity()
        
        import time
        #g = (time.time()/2.) % 1
        from math import sin, pi
        g = .5 + .5 * sin( 2*pi*( time.time()/2. ) )
        glClearColor( g,g,g,0 )
        glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT )
        
        '''
        glColor3f( 0,0,1 )
        glBegin( GL_QUADS )
        glVertex2f( -.5, -.5 )
        glVertex2f( .5, -.5 )
        glVertex2f( .5, .5 )
        glVertex2f( -.5, .5 )
        glEnd()
        '''
        
        self.SwapBuffers()
    
    def SaveToFile( self, fname, front_buffer = True ):
        '''
        Records the current framebuffer, front or back depending on the 'front_buffer' parameter,
        to the file named 'fname'.
        The image format is induced by the .extension of 'fname'.
        
        TODO: This function calls SetCurrent() on the glcanvas, and I don't know how to push/pop
              this bit of state.  Therefore, be careful when calling this inside of an opengl
              drawing function.
        '''
        
        self.SetCurrent( self.context )
        
        if front_buffer:
            read_buffer = glGetIntegerv( GL_READ_BUFFER )
            glReadBuffer( GL_FRONT )
        
        view = glGetIntegerv( GL_VIEWPORT )
        data = glReadPixels( view[0], view[1], view[2], view[3], GL_RGBA, GL_UNSIGNED_BYTE )
        
        if front_buffer: glReadBuffer( read_buffer )
        
        from PIL import Image
        img = Image.fromstring( 'RGBA', ( view[2], view[3] ), data )
        img = img.transpose( Image.FLIP_TOP_BOTTOM )
        img.save( fname )
    
    def event2screen( self, pos ):
        return ( pos[0], self.size[1] - 1 - pos[1] )
    
    def event2eye( self, pos ):
        return self.screen2eye( self.event2screen( pos ) )
    
    def event2world( self, pos ):
        return self.screen2world( self.event2screen( pos ) )
    
    def screen2eye( self, pos ):
        half_width = .5 * self.size[0]
        half_height = .5 * self.size[1]
        
        if self.size.width >= self.size.height:
            aspect = self.size.width / float( self.size.height )
            
            x = aspect * ( pos[0] - half_width ) / half_width
            y = ( pos[1] - half_height ) / half_height
        else:
            aspect = self.size.height / float( self.size.width )
            
            x = ( pos[0] - half_width ) / half_width
            y = aspect * ( pos[1] - half_height ) / half_height
        
        return x, y
    
    def screen2world( self, pos ):
        # print '===screen2world in:', pos
        #from helpers import solve_for_fixed_points
        #result = solve_for_fixed_points( [pos], ((0,0), self.size), ((-1,1), (1,-1)) )[0]
        result = (
            #-1. + 2 * float( pos[0] )/self.size[0],
            #-1. + 2 * float( pos[1] )/self.size[1]
            
            self.world_bounds[0] + (self.world_bounds[1]-self.world_bounds[0]) * float( pos[0] )/self.size[0],
            self.world_bounds[2] + (self.world_bounds[3]-self.world_bounds[2]) * float( pos[1] )/self.size[1]
            )
        #print '===screen2world(', pos, '):', result
        return result
    
    def world2screen( self, pos ):
        # print '===world2screen in:', pos
        result = (
            #self.size[0] * .5 * ( pos[0] + 1. ),
            #self.size[1] * .5 * ( pos[1] + 1. )
            
            self.size[0] * float( pos[0] - self.world_bounds[0] ) / (self.world_bounds[1]-self.world_bounds[0]),
            self.size[1] * float( pos[1] - self.world_bounds[2] ) / (self.world_bounds[3]-self.world_bounds[2])
            )
        #print '===world2screen(', pos, '):', result
        return result
    
    def glLineWidthWorld( self, world_width ):
        '''
        Given a line width in world-space, calls glLineWidth() with the appropriate pixel width.
        '''
        ## These scales should be the same, since the world_bounds and the pixel size of the window
        ## have the same aspect ratio.
        xscale = self.size[0] / float(self.world_bounds[1]-self.world_bounds[0])
        yscale = self.size[1] / float(self.world_bounds[3]-self.world_bounds[2])
        average_scale = .5 * ( xscale + yscale )
        glLineWidth( average_scale * world_width )


class ExampleFrame( wx.Frame ):
    def __init__(
            self, parent, id, title,
            pos = wx.DefaultPosition, size = wx.DefaultSize, style = wx.DEFAULT_FRAME_STYLE
            ):
        wx.Frame.__init__( self, parent, id, title, pos = pos, size = size, style = style )
        
        self.canvas = BasicGLCanvas( self, debugging = True )
        
        ## Forces a redisplay 30 times per second.
        from wx_util import Timers
        Timers.RequestTimer( lambda: self.canvas.Refresh(), 1000./30., False )
        
        '''
        menubar = wx.MenuBar()
        
        filem = wx.Menu()
        
        display = wx.Menu()
        
        item = display.Append( wx.ID_ANY, "Show Target" )
        def set_which( which ): self.canvas.draw_which = which
        self.Bind( wx.EVT_MENU, lambda evt: set_which( 'target' ), item )
        item = display.Append( wx.ID_ANY, "Show Undeformed" )
        self.Bind( wx.EVT_MENU, lambda evt: set_which( 'undeformed' ), item )
        item = display.Append( wx.ID_ANY, "Show Verify" )
        self.Bind( wx.EVT_MENU, lambda evt: set_which( 'verify' ), item )
        
        #display.AppendSeparator()
        #item = display.Append( wx.ID_ANY, "Clea&r User Points\tctrl-r" )
        #self.Bind( wx.EVT_MENU, self.OnClearPoints, item )
        
        modify = wx.Menu()
        
        item = modify.Append( wx.ID_ANY, "Increase Expansion\tctrl-=" )
        self.Bind( wx.EVT_MENU, lambda evt: (self.canvas.expansion_increase(), self.canvas.dough.update_solutions()), item )
        item = modify.Append( wx.ID_ANY, "Decrease Expansion\tctrl--" )
        self.Bind( wx.EVT_MENU, lambda evt: (self.canvas.expansion_decrease(), self.canvas.dough.update_solutions()), item )
        modify.AppendSeparator()
        item = modify.Append( wx.ID_ANY, "Toggle Gravity\tctrl-g" )
        self.Bind( wx.EVT_MENU, self.OnToggleGravity, item )
        modify.AppendSeparator()
        item = modify.Append( wx.ID_ANY, "Sketch Boundary\tctrl-d" )
        self.Bind( wx.EVT_MENU, self.OnMenuDrawBoundary, item )
        #modify.AppendSeparator()
        item = modify.Append( wx.ID_ANY, "Update Solutions\tctrl-u" )
        self.Bind( wx.EVT_MENU, lambda evt: self.canvas.dough.update_solutions(), item )
        
        
        menubar.Append( filem, "&File" )
        menubar.Append( display, "&Display" )
        menubar.Append( modify, "&Modify" )
        self.SetMenuBar( menubar )
        '''
        
        self.Show( True )


def load_PyCrust( app ):
    from wx.py.crust import CrustFrame
    crust = CrustFrame()
    crust.SetSize((750,525))
    crust.Show(True)
    crust.shell.interp.locals['app'] = app
    return crust

def main():
    app = wx.App()
    
    ## load the PyCrust shell
    ## NOTE: This generates the following errors and is very slow switching between the app window
    ##       and the PyCrust window.
    ## Sun Jun 15 20:19:36 yotronagon.local Python[32825] <Error>: CGContextRestoreGState: invalid context
    # crust = load_PyCrust( app )
    
    #frame = wx.Frame( None, -1, 'simple view', size=(512,512) )
    #canvas = BasicGLCanvas( frame )
    #frame.Show(True)
    
    frame = ExampleFrame( None, -1, 'example view', size = (512,512) )
    
    app.MainLoop()
    app.Destroy()

if __name__ == '__main__': main()
