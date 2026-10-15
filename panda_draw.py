import direct.directbase.DirectStart
from pandac.PandaModules import *
from direct.showbase.DirectObject import DirectObject
from direct.task.Task import Task

dbg_evts = False

class EventDispatcher( DirectObject ):
    def __init__( self, evt_handler = None ):
        self.running = True
        
        self.accept( 'mouse1', self.mouse_down )
        self.accept( 'mouse1-up', self.mouse_up )
        taskMgr.add( self.mouse_sampler, 'mouse-sample' )
        self.last_mouse_pos = None
        self.mouse_is_down = False
        
        self.event_handler = evt_handler
    
    def destroy( self ):
        self.ignoreAll()
        taskMgr.remove( 'mouse-sample' )
    
    def pos_is_ok_and_new( self, pos ):
        if pos is None: return False
        if self.last_mouse_pos is None: return True
        #print pos, self.last_mouse_pos
        if tuple( pos ) != tuple( self.last_mouse_pos ): return True
    
    def mouse_down( self ):
        if dbg_evts: print('mouse_down')
        self.mouse_is_down = True
        
        pos = get_mouse_pos()
        if pos is not None:
            if dbg_evts: print('pos ok')
            self.last_mouse_pos = pos
            if self.event_handler is not None:
                self.event_handler.doMouseDown( 0, pos )
    
    def mouse_up( self ):
        if dbg_evts: print('mouse_up')
        self.mouse_is_down = False
        
        pos = get_mouse_pos()
        if pos is not None:
            if dbg_evts: print('pos ok')
            self.last_mouse_pos = pos
            if self.event_handler is not None:
                self.event_handler.doMouseUp( 0, pos )
    
    def mouse_sampler( self, task ):
        pos = get_mouse_pos()
        if self.pos_is_ok_and_new( pos ):
            if dbg_evts: print('mouse_sampler dragging is', self.mouse_is_down)
            self.last_mouse_pos = pos
            if self.event_handler is not None:
                if self.mouse_is_down: self.event_handler.doMouseDrag( pos )
                else: self.event_handler.doMouseMove( pos )
        
        return Task.cont

def get_mouse_pos():
    if base.mouseWatcherNode.hasMouse():
        mpos = base.mouseWatcherNode.getMouse() #get the mouse position
        return ( mpos.getX(), mpos.getY() )
    else:
        return None

__mouse_does_rotate_camera = True
def get_mouse_rotates_camera():
    return __mouse_does_rotate_camera
def set_mouse_rotates_camera( whether ):
    if whether:
        mat=Mat4(camera.getMat())
        mat.invertInPlace()
        base.mouseInterfaceNode.setMat(mat)
        base.enableMouse()
    else:
        base.disableMouse()
    
    __mouse_does_rotate_camera = whether

class EventHandler( object ):
    def doMouseDown( self, button, pos ): pass
    def doMouseUp( self, button, pos ): pass
    def doMouseDrag( self, pos ): pass
    def doMouseMove( self, pos ): pass


class Sketcher( EventHandler ):
    def __init__( self ):
        self.n = 0
        self.curve = []
    
    def doMouseDown( self, button, pos ):
        if dbg_evts: print('sketcher.mdown')
        self.curve = [ pos ]
    
    def doMouseDrag( self, pos ):
        if dbg_evts: print('sketcher.mdrag')
        
        kDragThreshold = .01
        from numpy import asarray
        from myarray import mag
        if len( self.curve ) > 0 and mag( asarray( self.curve[-1] ) - asarray( pos ) ) < kDragThreshold: return
        print(pos)
        self.curve.append( pos )
    
    def doMouseUp( self, button, pos ):
        if dbg_evts: print('sketcher.mup')
        
        '''
        if len( self.curve ) < 3:
            set_mouse_rotates_camera( not get_mouse_rotates_camera() )
            print 'Click! Mouse rotates camera is now', get_mouse_rotates_camera()
            return
        '''
        
        import fibermesh
        import panda_util
        
        #curve = fibermesh.smooth_curve( self.curve )
        curve = self.curve
        print(curve)
        #mesh = fibermesh.inflate( curve )
        mesh = fibermesh.test2()
        node = panda_util.create_TrianglesNode_for_TriMesh( mesh, 'inflated_stroke_%s' % (self.n,) )
        NodePath( node ).reparentTo( render )
        
        self.n += 1
        self.curve = []
        
        set_mouse_rotates_camera( True )

def setMouseEnabledPos( position, rotation ):
    '''
    From http://panda3d.org/phpbb2/viewtopic.php?t=3867&highlight=geomvertexdata
    '''
    from pandac.PandaModules import Mat4
    # set a default camera position
    camera.setPosHpr( position, rotation )
    #camera.setPos(position)
    camera.lookAt(render)
    #camera.setHpr(camera, rotation)
    cameraMat=Mat4(camera.getMat())
    cameraMat.invertInPlace()
    base.mouseInterfaceNode.setMat(cameraMat) 

def init():
    set_mouse_rotates_camera( False )
    
    ## An orthographic camera to match the 2d camera
    lens = OrthographicLens()
    lens.setFilmSize( 2, 2 )
    lens.setNearFar( -1000, 1000 )
    base.cam.node().setLens( lens )
    
    # camera.setPosHpr( 0, 0, 25, 0, -90, 0 )  #Place the camera
    camera.setPosHpr( 0, 0, 25, 0, 0, 0 )  #Place the camera
    # setMouseEnabledPos( (0,0,25), (0, -90, 0) )
    # setMouseEnabledPos( (0,0,25), (0, 0, 0) )
    
    #sketcher = Sketcher()
    #e = EventDispatcher( sketcher )
    
def init_test():
    from mesh_gui import TriMesh
    import panda_util
    mesh = TriMesh.FromOBJ_FileName( 'optimized9.obj' )
    #print mesh.vs
    #print mesh.faces
    node = panda_util.create_TrianglesNode_for_TriMesh( mesh, 'optimized9' )
    NodePath( node ).reparentTo( render )
    
    camera.lookAt(render)
    camera.setPosHpr( 0, 0, 25, 0, -90, 0 )
    
    run()

def main():
    #init()
    init_test()

if __name__ == '__main__': main()
