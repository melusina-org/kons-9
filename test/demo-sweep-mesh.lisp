(in-package #:kons-9)

#|
A SWEEP-MESH is a subclass of UV-MESH and DEPENDENCY-NODE-MIXIN. It uses the
SWEEP-EXTRUDE-AUX function of UV-MESH to generate and update its mesh. The
updating is done using the DEPENDENCY-NODE-MIXIN class.

The demos below demonstrate examples of generating and updating SWEEP-MESH
instances.
|#

#|
(Demo 01 sweep-mesh) simple sweep-mesh =========================================

Create a SWEEP-MESH and modify its slots. The mesh automatically updates.
|#

(format t "  sweep-mesh 1...~%") (finish-output)

(with-clear-scene
  (defparameter *profile* (make-circle 0.8 4))
  (defparameter *path* (make-sine-curve 360 1 4 1 32))
  (defparameter *mesh* (make-sweep-mesh *profile* 0 *path* 0 :twist (* 2 pi) :taper 0.0))
  (add-shape *scene* *mesh*))

;;; modify slots and shape will change due to DEPENDENCY-NODE-MIXIN updating
(setf (num-segments *profile*) 8)
(setf (num-segments *path*) 64)
(setf (taper *mesh*) 1.0)

#|
(Demo 02 sweep-mesh) animate a sweep-mesh slot =================================

Create a SWEEP-MESH and modify animate a slot. The mesh automatically updates.
|#

(format t "  sweep-mesh 2...~%") (finish-output)

(with-clear-scene
  (defparameter *profile* (make-circle 1.2 4))
  (defparameter *path* (make-sine-curve 360 1 4 1 32))
  (defparameter *mesh* (make-sweep-mesh *profile* 0 *path* 0 :twist (* 2 pi) :taper 0.0))
  (let ((anim (make-instance 'animator
                             :setup-fn (lambda () (setf (num-segments *profile*) 4) nil)
                             :update-fn (lambda () (incf (num-segments *profile*))))))
    (add-motion *scene* anim)
    (add-shape *scene* *mesh*)))
;;; hold down space key in 3D view to run animation

#|
(Demo 03 sweep-mesh) animate a sweep-mesh-group ================================

Create a SWEEP-MESH-GROUP along a PARTICLE-SYSTEM. The mesh automatically grows
with the particle system.
|#

(format t "  sweep-mesh 3...~%") (finish-output)

(with-clear-scene
  (let* (;(p-gen (make-grid-uv-mesh 8 8 24 24))
         (p-sys (make-particle-system-from-point (p! 0 0 0) 10 (p! -.2 .4 -.2) (p! .2 .5 .2)
                                                 'particle
                                                 :update-angle (range-float 20.0 10.0)))
         (sweep-mesh-group (make-sweep-mesh-group (make-circle 0.5 6)
                                                  p-sys
                                                  :taper 0.0 :twist 2pi)))
    (setf (name sweep-mesh-group) 'sweep-group)
    (add-shape *scene* sweep-mesh-group)
    (add-motion *scene* p-sys)))
;;; hold down space key in 3D view to run animation

#|
(Demo 04 sweep-mesh) animate a sweep-mesh-group ================================

Create a SWEEP-MESH-GROUP along a PARTICLE-SYSTEM. The mesh automatically grows
with the particle system.
|#

(format t "  sweep-mesh 4...~%") (finish-output)

(with-clear-scene
  (let* ((p-source (make-cube-sphere 2 2))
         (p-sys (make-particle-system-from-point-source p-source
                                                        (lambda (v) (p:scale v 0.2))
                                                        'particle
                                      :life-span (round (rand2 5 10))
                                      :update-angle (range-float 20.0 10.0)))
         (sweep-mesh-group (make-sweep-mesh-group (make-circle 0.1 6) p-sys
                                                  :taper 1.0 :twist 0.0)))
    (add-shape *scene* p-source)
    (add-shape *scene* sweep-mesh-group)
    (add-motion *scene* p-sys)))
;;; hold down space key in 3D view to run animation

#|
(Demo 05 sweep-mesh) animate a sweep-mesh-group ================================

Create a SWEEP-MESH-GROUP along a PARTICLE-SYSTEM growing along a POLYHEDRON.
The mesh automatically grows with the particle system.
|#

(format t "  sweep-mesh 5...~%") (finish-output)

(with-clear-scene
  (let* ((shape (triangulate-polyhedron (make-cube-sphere 6.0 3)))
         ;(shape (triangulate-polyhedron (make-cube 6)))
         ;(shape (triangulate-polyhedron (make-grid-uv-mesh 6 6 1 1)))
         (p-sys (make-particle-system-from-point (p! 0 3.0 0) 1 (p! -.5 0 -.5) (p! .5 0 .5)
                                                 'climbing-particle
                                                 :support-polyh shape
                                                 :update-angle (range-float 20.0 10.0)
                                                 :life-span -1))
         (sweep-mesh-group (make-sweep-mesh-group (make-circle 0.2 6) p-sys
                                                  :taper 1.0 :twist 0.0)))
    (add-shape *scene* shape)
    (add-shape *scene* p-sys)
    (add-shape *scene* sweep-mesh-group)
    (add-motion *scene* p-sys)))
;;; hold down space key in 3D view to run animation -- gets slow, need to optimize
;;; suggestion: turn off filled display for a better view (TAB, D, 1)

(update-scene *scene* 240)               ;do update for batch testing

#|
(Demo 06 sweep-mesh) animate a sweep-mesh-group ================================

Create a SWEEP-MESH-GROUP along a PARTICLE-SYSTEM. The mesh automatically grows
with the particle system.
|#

(format t "  sweep-mesh 6...~%") (finish-output)

(with-clear-scene
  (let* ((p-sys (make-particle-system-from-point (p! 0 0 0) 10 (p! -.2 .4 -.2) (p! .2 .8 .2)
                                                 'dynamic-particle
                                                 :life-span -1 ;infinite life-span
                                                 :do-collisions? t
                                                 :force-fields (list (make-instance 'constant-force-field
                                                                                    :force-vector (p! 0 -0.05 0)))))
         (sweep-mesh-group (make-sweep-mesh-group (make-circle 0.5 6)
                                                  p-sys
                                                  :taper 0.0 :twist 2pi)))
    (setf (name sweep-mesh-group) 'sweep-group)
    (add-shape *scene* sweep-mesh-group)
    (add-motion *scene* p-sys)))
;;; hold down space key in 3D view to run animation

