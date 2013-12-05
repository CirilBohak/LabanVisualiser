var camera, scene, renderer;
var geometry, material, mesh, light;

init();
animate();

function init() {

    camera = new THREE.PerspectiveCamera( 75, window.innerWidth / window.innerHeight, 1, 10000 );
    camera.position.z = 1000;

    scene = new THREE.Scene();

    geometry = new THREE.SphereGeometry( 200, 50, 50 );
    //material = new THREE.MeshBasicMaterial( { color: 0x336699 } );
    material = new THREE.MeshPhongMaterial( { ambient: 0x020202, color: 0x336699, specular: 0x111111,  } );

    mesh = new THREE.Mesh( geometry, material );
    mesh.scale = new THREE.Vector3( 2.0, 2.0, 2.0 );
    scene.add( mesh );

    light = new THREE.PointLight(0xffffff, 2, 1000);
    light.position = new THREE.Vector3(200, 200, 1000);
    scene.add(light);

    renderer = new THREE.WebGLRenderer( { alpha: true, antialias: true } );
    renderer.setSize( window.innerWidth, window.innerHeight );

    document.body.appendChild( renderer.domElement );

}

function animate() {

    requestAnimationFrame( animate );

    mesh.rotation.x += 0.01;
    mesh.rotation.y += 0.02;

    renderer.render( scene, camera );

}