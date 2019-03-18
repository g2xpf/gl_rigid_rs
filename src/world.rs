use crate::circle::CircleBrander;
use crate::consts;
use crate::rectangle::RectangleBrander;
use crate::types::*;

pub struct World<'a> {
    pub circle_brander: CircleBrander<'a>,
    pub rectangle_brander: RectangleBrander<'a>,
    pub bodies: Vec<E>,
    pub uniform: GlobalUniform,

    pub collisions: Vec<(Collision)>,
}

impl<'a> World<'a> {
    pub fn new(
        circle_brander: CircleBrander<'a>,
        rectangle_brander: RectangleBrander<'a>,
        uniform: GlobalUniform,
    ) -> Self {
        World {
            circle_brander: circle_brander,
            rectangle_brander: rectangle_brander,
            bodies: Vec::new(),
            uniform: uniform,
            collisions: Vec::new(),
        }
    }

    pub fn brand<T>(&self, target: &mut T)
    where
        T: glium::Surface,
    {
        for body in &self.bodies {
            match body {
                E::Circle(circle) => {
                    if circle.visible == false {
                        continue;
                    }

                    target
                        .draw(
                            &self.circle_brander.vertex_buffer,
                            &self.circle_brander.index_buffer,
                            &self.circle_brander.program,
                            &uniform! {
                                r: circle.r,
                                angle: circle.angle,
                                radius: circle.radius,
                                scale: self.uniform.scale,
                                resolution: self.uniform.resolution,
                            },
                            &self.circle_brander.draw_parameter,
                        )
                        .unwrap();
                }
                E::Rectangle(rectangle) => {
                    if rectangle.visible == false {
                        continue;
                    }

                    target
                        .draw(
                            &self.rectangle_brander.vertex_buffer,
                            &self.rectangle_brander.index_buffer,
                            &self.rectangle_brander.program,
                            &uniform! {
                                r: rectangle.r,
                                wh: rectangle.wh,
                                scale: self.uniform.scale,
                            },
                            &self.circle_brander.draw_parameter,
                        )
                        .unwrap();
                }
            }
        }
    }

    pub fn step(&mut self) {
        for _ in 0..5 {
            self.init_constraints();
            self.detect_collisions();
            self.resolve_constraints();
        }

        self.integrate();
    }

    fn init_constraints(&mut self) {
        self.collisions = Vec::new();
    }

    fn detect_collisions(&mut self) {
        for (i, b1) in self.bodies.iter().enumerate() {
            for (j, b2) in self.bodies.iter().enumerate() {
                if i >= j {
                    continue;
                }
                match (&b1, &b2) {
                    (E::Circle(c1), E::Circle(c2)) => {
                        if c1.is_static && c2.is_static {
                            continue;
                        }
                        let dx = c2.r[0] - c1.r[0];
                        let dy = c2.r[1] - c1.r[1];
                        let dr1 = (c1.radius + c2.radius) * (c1.radius + c2.radius);
                        let dr2 = dx * dx + dy * dy;
                        let diff = dr1 - dr2;
                        if diff > 0.0 {
                            if dr2 == 0.0 {
                                continue;
                            }
                            let nx = -dx / dr2.sqrt();
                            let ny = -dy / dr2.sqrt();
                            self.collisions.push(Collision {
                                lhs: i,
                                rhs: j,
                                normal: [nx, ny],
                                depth: diff.sqrt(),
                                collision_point: [
                                    c1.r[0] - c1.radius * nx,
                                    c1.r[1] - c1.radius * ny,
                                ],
                            });
                        }
                    }
                    _ => (),
                }
            }
        }
    }

    fn resolve_constraints(&mut self) {
        for &Collision {
            lhs,
            rhs,
            normal,
            depth,
            collision_point,
        } in self.collisions.iter()
        {
            let force_n;
            let mut force_t;

            let nx = normal[0];
            let ny = normal[1];
            let tx = normal[1];
            let ty = -normal[0];

            let r1_cross_n;
            let r2_cross_n;
            let r1_cross_t;
            let r2_cross_t;
            match (&self.bodies[lhs], &self.bodies[rhs]) {
                (E::Circle(c1), E::Circle(c2)) => {
                    let r1x = collision_point[0] - c1.r[0];
                    let r1y = collision_point[1] - c1.r[1];
                    let r2x = collision_point[0] - c2.r[0];
                    let r2y = collision_point[1] - c2.r[1];

                    let rel_vx = (c2.v[0] - r2y * c2.omega) - (c1.v[0] - r1y * c1.omega);
                    let rel_vy = (c2.v[1] + r2x * c2.omega) - (c1.v[1] + r1x * c1.omega);

                    r1_cross_n = r1x * ny - r1y * nx;
                    r2_cross_n = r2x * ny - r2y * nx;
                    r1_cross_t = r1x * ty - r1y * tx;
                    r2_cross_t = r2x * ty - r2y * tx;

                    let rel_vn = rel_vx * nx + rel_vy * ny;
                    let rel_vt = rel_vx * tx + rel_vy * ty;

                    let mass_n = c1.inv_mass
                        + c2.inv_mass
                        + c1.inv_inertia * r1_cross_n * r1_cross_n
                        + c2.inv_inertia * r2_cross_n * r2_cross_n;
                    let mass_t = c1.inv_mass
                        + c2.inv_mass
                        + c1.inv_inertia * r1_cross_t * r1_cross_t
                        + c2.inv_inertia * r2_cross_t * r2_cross_t;

                    // calc normal force

                    // calc restitution
                    let e = if rel_vn < consts::RESTITUTION_BORDER {
                        0.0
                    } else {
                        consts::E
                    };
                    let mut rel_vn_prime = -e * rel_vn;
                    if rel_vn_prime > 0.0 {
                        rel_vn_prime = 0.0;
                    }

                    if depth >= consts::EQUILIBRIUM {
                        let separation_velocity =
                            -(depth - consts::EQUILIBRIUM) * consts::K / self.uniform.dt;
                        if rel_vn_prime > separation_velocity {
                            rel_vn_prime = separation_velocity;
                        }
                    }

                    force_n = -(rel_vn_prime - rel_vn) / mass_n;

                    // calc penalty force
                    // if depth >= EQUILIBRIUM {
                    //     let separation_force = K * (depth - EQUILIBRIUM) * self.uniform.dt;
                    //     if force_n < separation_force {
                    //         force_n = separation_force;
                    //     }
                    // }

                    force_t = -(0.0 - rel_vt) / mass_t;
                    let force_t_limit = force_n.abs() * consts::MU;
                    if force_t.abs() > force_t_limit {
                        force_t = force_t.signum() * force_t_limit;
                    }
                }
                _ => panic!("invalid collision"),
            }

            match &mut self.bodies[lhs] {
                E::Circle(c) => {
                    c.v[0] += force_n * nx * c.inv_mass;
                    c.v[1] += force_n * ny * c.inv_mass;
                    c.omega += force_n * r1_cross_n * c.inv_mass;
                    c.v[0] += force_t * tx * c.inv_mass;
                    c.v[1] += force_t * ty * c.inv_mass;
                    c.omega += force_t * r1_cross_t * c.inv_mass;
                }
                _ => (),
            }

            match &mut self.bodies[rhs] {
                E::Circle(c) => {
                    c.v[0] -= force_n * nx * c.inv_mass;
                    c.v[1] -= force_n * ny * c.inv_mass;
                    c.omega -= force_n * r2_cross_n * c.inv_mass;
                    c.v[0] -= force_t * tx * c.inv_mass;
                    c.v[1] -= force_t * ty * c.inv_mass;
                    c.omega -= force_t * r2_cross_t * c.inv_mass;
                }
                _ => (),
            }
        }
    }

    fn integrate(&mut self) {
        for body in &mut self.bodies {
            match body {
                E::Circle(circle) => {
                    if circle.is_static {
                        circle.a[0] = 0.0;
                        circle.a[1] = 0.0;
                        circle.v[0] = 0.0;
                        circle.v[1] = 0.0;
                        continue;
                    };
                    circle.a[0] = self.uniform.gravity[0];
                    circle.a[1] = self.uniform.gravity[1];
                    circle.v[0] += circle.a[0] * self.uniform.dt;
                    circle.v[1] += circle.a[1] * self.uniform.dt;
                    circle.r[0] += circle.v[0] * self.uniform.dt;
                    circle.r[1] += circle.v[1] * self.uniform.dt;
                    circle.angle += circle.omega * self.uniform.dt;
                }
                E::Rectangle(rectangle) => {
                    if rectangle.is_static {
                        rectangle.a[0] = 0.0;
                        rectangle.a[1] = 0.0;
                        rectangle.v[0] = 0.0;
                        rectangle.v[1] = 0.0;
                        continue;
                    };
                    rectangle.a[0] = self.uniform.gravity[0];
                    rectangle.a[1] = self.uniform.gravity[1];
                    rectangle.v[0] += rectangle.a[0] * self.uniform.dt;
                    rectangle.v[1] += rectangle.a[1] * self.uniform.dt;
                    rectangle.r[0] += rectangle.v[0] * self.uniform.dt;
                    rectangle.r[1] += rectangle.v[1] * self.uniform.dt;
                    rectangle.angle += rectangle.omega * self.uniform.dt;
                }
            }
        }
    }
}
