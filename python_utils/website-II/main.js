import * as THREE from 'three';
import { PLYLoader } from './PLYLoader.js';

import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { DragControls } from 'three/addons/controls/DragControls.js';

let container = document.getElementById("CanvasContainer");
let slider = document.getElementById("frameselector");

let stats, orbitControls, dragControls;

let camera, cameraTarget, scene, renderer;
let bubbles = new THREE.Group();
//bubbles.visible = false;

// make only the 
slider.oninput = (event) => {
    set_frame(event.target.value);
    document.getElementById('slidenumber').textContent = event.target.value;
    draw();
};


init();
//set_frame(10); // doesn't work right now, since ply's not loaded at this point
//animate();
draw();



let index_list = [];

function thingsToDo() {
    console.log('I do them now.');
    bubbles.children.forEach((child)=> {
        index_list.push(child.frame_index);
    });
    slider.min = Math.min.apply(Math,index_list);
    slider.max = Math.max.apply(Math,index_list);
    slider.value = 0;
    set_frame(0);
    draw();
}

const loader = new PLYLoader();

let fileInput = document.getElementById('ply-input');
let fileList = [];
let remaining = {};

fileInput.addEventListener('change', (e)=> { 
    
    fileList = [];
    remaining.val = fileInput.files.length;
    for(let i = 0; i < fileInput.files.length; i++) {
        fileList.push(fileInput.files[i]);
    }
    for(let i=0;i<fileList.length;i++) {

        let re = /(?:\.([^.]+))?$/;

        let filename = fileList[i].name;

        let ext = re.exec(filename)[1];

        switch(ext) {
            case 'ply':
                readPlyFile(loader,fileList[i],bubbles,()=>{
                    checkFinish(remaining,thingsToDo);
                });
            break;

            case 'csv':
                if(filename == 'times.csv') {
                    console.log('found the times!',fileList[i]);
                }
                checkFinish(remaining,thingsToDo);
            break;

            case undefined: // no break!
                console.log('unknown file type encountered');

            default:
                checkFinish(remaining,thingsToDo);
        }
        
        //console.log(fileList[i].name);
        //const url = window.URL.createObjectURL(fileList[i]);
        //console.log(url);
        //load_ply_geometry(loader,url,bubbles);

        
        
    }
    

});

function checkFinish(counter,callback) {
    if(counter.val <= 0) {
        console.log('error: counter below or equal zero.');
    } else {
        counter.val--;
        if(counter.val == 0) {
            callback();
        }
    }   
}

function readPlyFile(plyloader,file,group,callback) {
    const reader = new FileReader();
    reader.addEventListener('load', async function(event) {
        const contents = event.target.result;
        const geometry = plyloader.parse(contents);
        geometry.computeVertexNormals();

        if(geometry.attributes.color === undefined && 
            geometry.attributes.potential !== undefined) compute_colors(geometry,[-1,1],greenblue);

        const material = new THREE.MeshPhongMaterial( {
            color: 0xffffff,
            flatShading: true,
            vertexColors: true,
            shininess: 0,
            polygonOffset: true, // these three lines create some sort of offset, to avoid z-fighting
            polygonOffsetFactor: 1, // of color surface and wireframe
            polygonOffsetUnits: 1
        } );

        const wireframeMaterial = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: true, transparent: true } );

        let mesh = new THREE.Mesh( geometry, material );
        let wireframe = new THREE.Mesh( geometry, wireframeMaterial );
        mesh.add(wireframe);

        //mesh.castShadow = true;
        //mesh.receiveShadow = true;

        const rx = /-(.*).ply/g;
        const i_a = parseInt(rx.exec(file.name)[1]);
        mesh.frame_index = i_a;

        group.add(mesh);
        callback(mesh);
    });
    reader.readAsArrayBuffer(file);
}


orbitControls.addEventListener('change', draw);
window.addEventListener('resize', draw);



function init() {

    scene = new THREE.Scene(); 
    console.log(scene);
    scene.background = new THREE.Color( 0x777777 );
    scene.fog = new THREE.Fog( 0x777777, 10, 25 );
    let axes = new THREE.AxesHelper(1.5);
    axes.material = new THREE.LineBasicMaterial( {
        vertexColors: true,
        linewidth: 3
    } );
    scene.add(axes);
    scene.rotation.x = - Math.PI / 2;

    camera = new THREE.PerspectiveCamera( 30, container.clientWidth / container.clientHeight, 0.1, 1000 );

    camera.position.set(0,5,0);

    cameraTarget = new THREE.Vector3( 0, 0, 0 );


    // Ground

    const plane = new THREE.Mesh(
        new THREE.PlaneGeometry( 1000, 1000 ),
        new THREE.MeshPhongMaterial( { color: 0x999999, specular: 0x101010 } )
    );
    plane.position.z = - 5;
    scene.add( plane );

    plane.receiveShadow = true;


    

    // PLY file
    scene.add(bubbles);

    // Lights

    //scene.add( new THREE.HemisphereLight( 0x443333, 0x111122 ) );
    scene.add(new THREE.AmbientLight(0xffffff));

    //addShadowedLight( 1, 0, 0, 0xffffff, 1.35 );
    //addShadowedLight( 0.5, 1, - 1, 0xffaa00, 1 );

    // renderer

    renderer = new THREE.WebGLRenderer( { antialias: true } );
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( container.clientWidth, container.clientHeight );

    renderer.shadowMap.enabled = true;

    container.appendChild( renderer.domElement );

    // stats

    //stats = new Stats();
    //container.appendChild( stats.dom );

    // resize

    window.addEventListener( 'resize', onWindowResize );

    // controls

    orbitControls = new OrbitControls(camera, renderer.domElement);
    dragControls = new DragControls(scene,camera, renderer.domElement);

    dragControls.addEventListener('dragstart', function (event) {
        orbitControls.enabled = false
        event.object.material.opacity = 0.33
    });
    dragControls.addEventListener('dragend', function (event) {
        orbitControls.enabled = true
        event.object.material.opacity = 1
    });
}

function addShadowedLight( x, y, z, color, intensity ) {

    const directionalLight = new THREE.DirectionalLight( color, intensity );
    directionalLight.position.set( x, y, z );
    scene.add( directionalLight );

    directionalLight.castShadow = true;

    const d = 1;
    directionalLight.shadow.camera.left = - d;
    directionalLight.shadow.camera.right = d;
    directionalLight.shadow.camera.top = d;
    directionalLight.shadow.camera.bottom = - d;

    directionalLight.shadow.camera.near = 1;
    directionalLight.shadow.camera.far = 4;

    directionalLight.shadow.mapSize.width = 1024;
    directionalLight.shadow.mapSize.height = 1024;

    directionalLight.shadow.bias = - 0.001;

}

function onWindowResize() {

    camera.aspect = container.clientWidth / container.clientHeight;
    camera.updateProjectionMatrix();

    renderer.setSize( container.clientWidth, container.clientHeight );

}

function draw() {
    requestAnimationFrame(render);
}

function animate() {

    render();
    stats.update();
    requestAnimationFrame( animate );

}

function render() {

    //const timer = Date.now() * 0.0005;

    //camera.position.x = Math.sin( timer ) * 2.5;
    //camera.position.z = Math.cos( timer ) * 2.5;

    //camera.lookAt( cameraTarget );

    renderer.render( scene, camera );

}

function load_ply_geometry(loader,fname, group) {
    loader.load( fname, function ( geometry) {

        console.log(geometry);

        geometry.computeVertexNormals();

        if(geometry.attributes.color === undefined && 
            geometry.attributes.potential !== undefined) compute_colors(geometry,[-1,1],greenblue);

        const material = new THREE.MeshPhongMaterial( {
            color: 0xffffff,
            flatShading: true,
            vertexColors: true,
            shininess: 0,
            polygonOffset: true, // these three lines create some sort of offset, to avoid z-fighting
            polygonOffsetFactor: 1, // of color surface and wireframe
            polygonOffsetUnits: 1
        } );

        const wireframeMaterial = new THREE.MeshBasicMaterial( { color: 0x000000, wireframe: true, transparent: true } );

        let mesh = new THREE.Mesh( geometry, material );
        let wireframe = new THREE.Mesh( geometry, wireframeMaterial );
        mesh.add(wireframe);

        mesh.castShadow = true;
        mesh.receiveShadow = true;

        const rx = /-(.*).ply/g;
        const i_a = parseInt(rx.exec(fname)[1]);
        mesh.frame_index = i_a;

        group.add( mesh );

    },()=>{},(err)=>{
        console.log(err);
    });

}

function compute_colors(geometry,range,colormap) {
    let potential = geometry.attributes.potential;

    let colorbuffer = [];
    for(let i=0;i<potential.count;i++){
        let color = generate_color(potential.getX(i),range,colormap);  
        colorbuffer.push(color.r,color.g,color.b); 
    }

    geometry.setAttribute('color',new THREE.Float32BufferAttribute(colorbuffer, 3));
}

function generate_color(value,range,colormap) {
    return colormap((value-range[0])/(range[1]-range[0]));
}

function greenblue(t) {
    let color = {};
    color.r = 0;
    color.g = Math.max(Math.min(1-t,1),0);
    color.b = Math.max(Math.min(t,1),0);
    return color;
}

function set_frame(index) {
    bubbles.children.forEach(child => {
        if(child.frame_index == index)
            child.visible = true;
        else
            child.visible = false;
    });
}










let coldown = '#00ff00';
let colup   = '#0000ff';

//coldown = '#0c0887';
//colup = '#f0f921';

function get_real_value(value,boundaryA,boundaryB) {
    value = boundaryA + (boundaryB-boundaryA)*value;

    return value;
}


let colcont = document.getElementById("ColorControl");

let colorgradient = document.createElement("div");
colorgradient.id = "colorgradient";
colcont.appendChild(colorgradient);


let slider_res = 1000.0;

let sliderA = document.createElement("input");
sliderA.className = 'colorslider';
sliderA.min = 0;
sliderA.max = slider_res;
sliderA.value = 0;
sliderA.type = 'range';
colcont.appendChild(sliderA);

let sliderB = document.createElement("input");
sliderB.className = 'colorslider';
sliderB.min = 0;
sliderB.max = slider_res;
sliderB.value = slider_res;
sliderB.type = 'range';
colcont.appendChild(sliderB);


let inputA = document.createElement('input');
inputA.className = 'boundaryvalue';
inputA.value = -1;
inputA.step = 0.01;
inputA.type = 'number';
colcont.appendChild(inputA);

let inputB = document.createElement('input');
inputB.className = 'boundaryvalue';
inputB.value = 1;
inputB.step = 0.01;
inputB.type = 'number';
inputB.style.float = 'right';
colcont.appendChild(inputB);


let sliders = colcont.getElementsByTagName('input');

for (let elm of sliders) {
    elm.oninput = (event) => {

        let valA = sliderA.value/slider_res*100.0;
        let valB = sliderB.value/slider_res*100.0;

        if(valA < valB) {
            colorgradient.style.background = 'linear-gradient(to right, '+coldown+' '+valA+'%, '+colup+' '+valB+'%)';
        }else{
            colorgradient.style.background = 'linear-gradient(to right, '+colup+' '+valB+'%, '+coldown+' '+valA+'%)';
        }
    
        
    
        bubbles.children.forEach( bubble => {
            
            let geometry = bubble.geometry;

            if(geometry.attributes.potential !== undefined) {
                let range = [get_real_value(sliderA.value/slider_res,Number(inputA.value),Number(inputB.value)),
                                get_real_value(sliderB.value/slider_res,Number(inputA.value),Number(inputB.value))];
                
                
                compute_colors(geometry,range,greenblue);
            }

            
        });

        draw();
    }
}
            