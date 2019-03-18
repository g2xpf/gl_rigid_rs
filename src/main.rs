extern crate gl_rigid_rs;
extern crate glium;
extern crate rand;

fn main() {
    use glium::{glutin, Surface};

    let mut events_loop = glutin::EventsLoop::new();
    let window = glutin::WindowBuilder::new()
        .with_title("GLWindow")
        .with_dimensions(glutin::dpi::LogicalSize::new(640.0, 640.0));
    let context = glutin::ContextBuilder::new().with_vsync(true);
    let display = glium::Display::new(window, context, &events_loop).unwrap();

    let circle_brander = gl_rigid_rs::CircleBrander::new(
        &display,
        glium::DrawParameters {
            blend: glium::Blend {
                alpha: glium::BlendingFunction::Addition {
                    source: glium::LinearBlendingFactor::SourceAlpha,
                    destination: glium::LinearBlendingFactor::OneMinusSourceAlpha,
                },
                ..Default::default()
            },
            ..Default::default()
        },
    );
    let rectangle_brander = gl_rigid_rs::RectangleBrander::new(&display, Default::default());
    let mut world = gl_rigid_rs::World::new(
        circle_brander,
        rectangle_brander,
        gl_rigid_rs::GlobalUniform {
            dt: 0.016,
            gravity: [0.0, -9.8],
            resolution: [640f32, 640f32],
            scale: 0.05,
        },
    );

    for i in 0..100 {
        let x = (i % 3) as f32;
        let x = 2.0 * x - 1.0;
        let y = (i / 3) as f32;
        let y = 2.0 * y - 1.0;
        world
            .bodies
            .push(gl_rigid_rs::E::Circle(gl_rigid_rs::Circle::new(
                0.6,
                [x, y + 10.0],
                [0.0, 0.0],
                false,
            )));
    }
    world
        .bodies
        .push(gl_rigid_rs::E::Circle(gl_rigid_rs::Circle::new(
            6.0,
            [10.0, -9.0],
            [0.0, 0.0],
            true,
        )));
    world
        .bodies
        .push(gl_rigid_rs::E::Circle(gl_rigid_rs::Circle::new(
            6.0,
            [-10.0, -9.0],
            [0.0, 0.0],
            true,
        )));
    world
        .bodies
        .push(gl_rigid_rs::E::Circle(gl_rigid_rs::Circle::new(
            10.0,
            [0.0, -20.0],
            [0.0, 0.0],
            true,
        )));

    let mut window_should_closed = false;
    while !window_should_closed {
        let mut target = display.draw();
        target.clear_color(0.03, 0.03, 0.03, 0.4);
        world.brand(&mut target);
        target.finish().unwrap();

        world.step();
        events_loop.poll_events(|ev| match ev {
            glutin::Event::WindowEvent { event, .. } => match event {
                glutin::WindowEvent::CloseRequested
                | glutin::WindowEvent::KeyboardInput {
                    input:
                        glutin::KeyboardInput {
                            scancode: 16,
                            state: glutin::ElementState::Pressed,
                            ..
                        },
                    ..
                } => window_should_closed = true,
                _ => (),
            },
            _ => (),
        });
    }
}
