
var canvas;
//var gl;
var crosssectionscanvas;
var ctx;
var voxdata;


var basicCubeData = {
	vertices : [
	/*
		0,0,0,	//0
		0,0,1,	//1
		0,1,0,	//2
		0,1,1,	//3
		1,0,0,	//4
		1,0,1,	//5
		1,1,0,	//6
		1,1,1],	//7
		*/
		
		0,0,0,	//0
		0,0,2,	//1
		0,2,0,	//2
		0,2,2,	//3
		2,0,0,	//4
		2,0,2,	//5
		2,2,0,	//6
		2,2,2],	//7
		
	indices: [	//these faces "inside-out" for easier testing
		0,2,1,	//z=0
		3,1,2,
				
		4,5,6,	//z=-1
		5,7,6,
		
		0,4,6,	//x=0
		0,6,2,
		
		1,3,7,	//x=1
		1,7,5,	
		
		0,1,4,	//y=0
		4,1,5,
		
		2,6,3,	//y=1
		6,7,3,
	]
}



var basicCubeBuffers={};
var voxData={};
var voxBuffers={};

function initBuffers(){

	loadBufferData(basicCubeBuffers, basicCubeData);
	
	voxBuffers["stupid"]={};
	voxBuffers["sparse"]={};
	voxBuffers["sparseWithNormals"]={};
	loadBufferData(voxBuffers["stupid"], voxData["stupid"]);
		
	loadBufferData(voxBuffers["sparse"], { vertices:voxData["sparse"].vertices, indices:voxData["sparse"].indices});	//copy all but normals!
	loadBufferData(voxBuffers["sparseWithNormals"], voxData["sparse"]);

	function bufferArrayData(buffer, arr, size){
		gl.bindBuffer(gl.ARRAY_BUFFER, buffer);
		gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(arr), gl.STATIC_DRAW);
		buffer.itemSize = size;
		buffer.numItems = arr.length / size;
		console.log("buffered. numitems: " + buffer.numItems);
	}
	
	function loadBufferData(bufferObj, sourceData){
		bufferObj.vertexPositionBuffer = gl.createBuffer();
		bufferArrayData(bufferObj.vertexPositionBuffer, sourceData.vertices, 3);
		//stuff about normals etc present in 3-sph project got this from, removed here. 
		
		if (sourceData.normals){
			bufferObj.vertexNormalBuffer = gl.createBuffer();
			bufferArrayData(bufferObj.vertexNormalBuffer, sourceData.normals, 3);
		}
		if (sourceData.colors){
			bufferObj.vertexColorBuffer = gl.createBuffer();
			bufferArrayData(bufferObj.vertexColorBuffer, sourceData.colors, 3);
		}
		
		//triangles rather than strip, but no big deal- frag shader does most of the work!
		bufferObj.vertexIndexBuffer = gl.createBuffer();
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bufferObj.vertexIndexBuffer);
		gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(sourceData.indices), gl.STATIC_DRAW);
		bufferObj.vertexIndexBuffer.itemSize = 3;
		bufferObj.vertexIndexBuffer.numItems = sourceData.indices.length;
	}
}


var shaderPrograms={};
function initShaders(){
	shaderPrograms.simple = loadShader( "shader-simple-vs", "shader-simple-fs",{
					attributes:["aVertexPosition"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
					
	shaderPrograms.withNormals = loadShader( "shader-withnormals-vs", "shader-color-fs",{
					attributes:["aVertexPosition","aVertexNormal"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
	
	shaderPrograms.withNormalsAndColor = loadShader( "shader-withnormals-andcolor-vs", "shader-color-fs",{
					attributes:["aVertexPosition","aVertexNormal","aVertexColor"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
}

var mvMatrix = mat4.create();
//mat4.identity(mvMatrix);
var pMatrix = mat4.create();
mat4.identity(pMatrix);
//mat4.translate(mvMatrix,vec3.fromArray([0,0,-10])); //glmatix expects own vector type


var playerPosition = [0,-0.5,-3];	//should be contained within playerMatrix. TODO sort this (code here consistent with webglwideanglecamera project)

mvMatrix[14]=-3;	//move back to look at thing (is this camera or thing position?)
mvMatrix[13]=-0.5;
mvMatrix[12]=0;

var playerMatrix = mat4.identity();

mat4.rotateX(playerMatrix, -0.5);	//rads
//mat4.rotateZ(playerMatrix, 0.2);
rotatePlayer([0,0,0]);	//hack to cause position update in matrix?



//mat4.set(mvMatrix, playerMatrix);	//??

var mainCamFov = 85;	//degrees.

var mouseInfo = {
	x:0,
	y:0,
	dragging: false,
	lastPointingDir:{},
	currentPointingDir:{x:0,y:0,z:1,w:1}
};


voxdata = [];

var guiParams={
	box:true,
	stupid:false,
	sparse:false,
	sparseWithNormals:true,	//TODO with/without color
	constrainToBox:false
};

var stats;



//voxel data. maybe better to have octree or similar. for now, use array of arrays. implementation maybe sparse. maybe good to use bits in a number (32 bit?)
var blocksize = 64;	//64x64x64 = 512x512, which seems manageable to draw slices onto canvas
var canvassize = 512;

function loadColData(data, object){
	//prep collision data for teapot
	var faces = [];
	var indices = object.indices;
	for (var ii=0;ii<indices.length;ii+=3){
		faces.push([indices[ii],indices[ii+1],indices[ii+2]]);
	}
	data.faces=faces;
	
	var tpVerts = object.vertices;
	var verts=[];
	for (var ii=0;ii<tpVerts.length;ii+=3){
		verts.push([tpVerts[ii],tpVerts[ii+1],tpVerts[ii+2]]);
	}
	data.verts=verts;
}

function loadBlenderExport(meshToLoad){
	return {
		vertices: meshToLoad.vertices,
		//normals: meshToLoad.normals,
		//uvcoords: meshToLoad.texturecoords?meshToLoad.texturecoords[0]:false,
		indices: [].concat.apply([],meshToLoad.faces)	//trick from https://www.youtube.com/watch?v=sM9n73-HiNA t~ 28:30
	}	
};


// do in n^2 checks if cast a series of parallel lines, assuming each starting from outside object. sort collisions, step through cells
var startGenFuncTime = Date.now();
var fromPolyModelFunctionFast = (function generateFromPolyModelFunctionFast(){
	
	var teapotColData={};	
	var teapotObject = loadBlenderExport(teapotData);	//isn't actually a blender export - just a obj json
	loadColData(teapotColData, teapotObject);
	
	var faces = teapotColData.faces;
	var verts = teapotColData.verts;
	
	//precalculate face normal data
	var facePlaneData=[];
	for (var ff=0;ff<faces.length;ff++){
		var thisTri = faces[ff];
		
		var triPoints = [verts[thisTri[0]], verts[thisTri[1]], verts[thisTri[2]]];
		var edgeVecs = [[triPoints[1][0] - triPoints[0][0], triPoints[1][1] - triPoints[0][1], triPoints[1][2] - triPoints[0][2]],
						[triPoints[2][0] - triPoints[0][0], triPoints[2][1] - triPoints[0][1], triPoints[2][2] - triPoints[0][2]]];
		var normVector = [edgeVecs[1][1]*edgeVecs[0][2] - edgeVecs[1][2]*edgeVecs[0][1],
							edgeVecs[1][2]*edgeVecs[0][0] - edgeVecs[1][0]*edgeVecs[0][2],
							edgeVecs[1][0]*edgeVecs[0][1] - edgeVecs[1][1]*edgeVecs[0][0]];
		
		facePlaneData.push({centreZ:triPoints[0][2] + ((normVector[0]*triPoints[0][0] + normVector[1]*triPoints[0][1])/normVector[2]),
							gradX: -normVector[0]/normVector[2],	
							gradY: -normVector[1]/normVector[2]});	//TODO check that gradients aren't /0
	}
	
	//pregenerate and store data.
	//to use some assumptions used to generate smooth data, return bilinear smoothed data
	var myVoxData = [];
	for (var ii=0;ii<blocksize;ii++){
		var slicedata = [];
		myVoxData.push(slicedata);
		for (var jj=0;jj<blocksize;jj++){
			var stripdata = [];
			slicedata.push(stripdata);
			var collisionData = checkForCollisions([(ii-30)/31,(jj-16)/31,-1],[(ii-30)/31,(jj-16)/31,1]);	//TODO rotate teapot 45 deg to fit into box better
			collisionData.sort(function(a,b){return a.z<b.z;});	//sort collision data.
			var zidx=0;
			var nextCollision;
			var fill=-1;
			while (nextCollision = collisionData.pop()){
				while (zidx<nextCollision.z){
					stripdata[zidx++]=fill;
				}
				fill = nextCollision.fill;	//transition to filled ie face downward or unfilled ie face up
			}
			while (zidx<blocksize){
				stripdata[zidx++]=fill;
			}
		}
	}
	console.log("myVoxData:");
	console.log(myVoxData);
	
	function checkForCollisions(raystart, rayend){
		var collisionPoints = [];
		var rayVec=[ rayend[0]-raystart[0], rayend[1]-raystart[1], rayend[2]-raystart[2] ];
	
		for (var ff=0;ff<faces.length;ff++){
			var thisTri = faces[ff];
			
			var signsSum=0;
			for (var vv=0;vv<3;vv++){
				var thisVert = verts[thisTri[vv]];
				var nextv = (vv+1) % 3;
				var nextVert = verts[thisTri[nextv]];
				
				/*
				var displacement = [raystart[0]-thisVert[0],raystart[1]-thisVert[1],raystart[2]-thisVert[2]];
				var edgevector = [nextVert[0]-thisVert[0],nextVert[1]-thisVert[1],nextVert[2]-thisVert[2]];
				var crossProd = [ displacement[1]*edgevector[2] - displacement[2]*edgevector[1],
									displacement[2]*edgevector[0] - displacement[0]*edgevector[2],
									displacement[0]*edgevector[1] - displacement[1]*edgevector[0]];
				var dotProd = crossProd[0]*rayVec[0] + crossProd[1]*rayVec[1] + crossProd[2]*rayVec[2];
				*/
				
				//simplify above given that rayVec x,y components =0
				var displacement = [raystart[0]-thisVert[0],raystart[1]-thisVert[1]];
				var edgevector = [nextVert[0]-thisVert[0],nextVert[1]-thisVert[1]];
				var crossProdZ = displacement[0]*edgevector[1] - displacement[1]*edgevector[0];
				var dotProd = crossProdZ*rayVec[2];
				
				//signsSum+=dotProd/Math.abs(dotProd);
				
				//use vertex order for handedness. might still have issues with ray through vertex.
				//TODO make this more efficient (takes ~ 20% longer than simple signsum calculation without if.
				if (thisTri[vv] > thisTri[nextv]){
					signsSum+=dotProd>=0 ? 1 : -1;
				}else{
					signsSum+=dotProd>0 ? 1 : -1;
				}
			}
			
			if (Math.abs(signsSum)>2.5){	//should need 3, but unsure how numerical error works
				//point on plane P
				//plane normal N 
				//looking for z for a line of fixed x,y
				// equation of plane q-p dot N = 0
				// -> qz-pz = (Nx*(qx-px) + Ny*(qy-py) ) / Nz
				/*
				var triPoints = [verts[thisTri[0]], verts[thisTri[1]], verts[thisTri[2]]];
				var edgeVecs = [[triPoints[1][0] - triPoints[0][0], triPoints[1][1] - triPoints[0][1], triPoints[1][2] - triPoints[0][2]],
								[triPoints[2][0] - triPoints[0][0], triPoints[2][1] - triPoints[0][1], triPoints[2][2] - triPoints[0][2]]];
				var normVector = [edgeVecs[1][1]*edgeVecs[0][2] - edgeVecs[1][2]*edgeVecs[0][1],
									edgeVecs[1][2]*edgeVecs[0][0] - edgeVecs[1][0]*edgeVecs[0][2],
									edgeVecs[1][0]*edgeVecs[0][1] - edgeVecs[1][1]*edgeVecs[0][0]];
				var zCollisionPoint= triPoints[0][2] - ((normVector[0]*(raystart[0]-triPoints[0][0]) + normVector[1]*(raystart[1]-triPoints[0][1]))/normVector[2]);
						// by calculation, - should be + above, presumably some mistake. 
				*/
				//precalc all this stuff once for each tri of object (rather than inside here every time collides with a ray). basically want x,y gradient and height at x,y=0 for each (non vertical) face.
				
				var facePlaneDataPoint = facePlaneData[ff];
				var zCollisionPoint= facePlaneDataPoint.centreZ + facePlaneDataPoint.gradX*raystart[0] + facePlaneDataPoint.gradY*raystart[1];
				
				collisionPoints.push({z:(zCollisionPoint+1)*32, fill:signsSum<0?-1:1});	//TODO check polarity
			}
		}
		
		//return some dummy data
		//return [{fill:1,z:5},{fill:0,z:10}];
		return collisionPoints;
	}
	
	return function(ii,jj,kk){
		//ii+=0.5;	//??
		//jj+=0.5;	//??
		//kk+=0.5;	//??
		
		iif=Math.floor(ii);
		jjf=Math.floor(jj);
		kkf=Math.floor(kk);
		//return myVoxData[iif][jjf][kkf];	//todo bilinear filter
		
		//remainder
		var iir = ii-iif;
		var jjr = jj-jjf;
		var kkr = kk-kkf;
		
		//bilinear filter version (bit slower. could code dedicated gradient extraction more efficiently)
		if (iif<63 && jjf<63 & kkf<63){
			var sum=0;
			sum+= myVoxData[iif][jjf][kkf]*(1-iir)*(1-jjr)*(1-kkr);
			sum+= myVoxData[iif][jjf][kkf+1]*(1-iir)*(1-jjr)*kkr;
			sum+= myVoxData[iif][jjf+1][kkf]*(1-iir)*jjr*(1-kkr);
			sum+= myVoxData[iif][jjf+1][kkf+1]*(1-iir)*jjr*kkr;
			sum+= myVoxData[iif+1][jjf][kkf]*iir*(1-jjr)*(1-kkr);
			sum+= myVoxData[iif+1][jjf][kkf+1]*iir*(1-jjr)*kkr;
			sum+= myVoxData[iif+1][jjf+1][kkf]*iir*jjr*(1-kkr);
			sum+= myVoxData[iif+1][jjf+1][kkf+1]*iir*jjr*kkr;
			return sum;
		}else{
			//strictly there are more edge cases but don't care for now
			//maybe this never gets hit... 
			return myVoxData[63][63][63];
			console.log("EDGE CASE HIT");
		}
		
	}
})();
console.log("generated poly to vox func in " + (Date.now() - startGenFuncTime) + " ms.");

//code borrowed from collision-detect-test project
//this seems to be extremely inefficient - n^3 checks. 
var fromPolyModelFunction = (function generateFromPolyModelFunction(){
	var teapotColData={};	
	var teapotObject = loadBlenderExport(teapotData);	//isn't actually a blender export - just a obj json
	loadColData(teapotColData, teapotObject);
	
	function checkForCollisions(raystart, rayend){
		var collisionPoints = [];
		//cross product rayVec with edge to ray, dot with edge vec
		//this is triple product, so can do in different order. see http://www.mathphysics.com/pde/vectorid.html
		
		var rayVec=[ rayend[0]-raystart[0], rayend[1]-raystart[1], rayend[2]-raystart[2] ];
	
		var faces = teapotColData.faces;
		var verts = teapotColData.verts;
		
		collidedFaces=0;
		
		for (var ff=0;ff<faces.length;ff++){
			var thisTri = faces[ff];
			
			//var isInside = true;
			var signsSum=0;
			for (var vv=0;vv<3;vv++){
				var thisVert = verts[thisTri[vv]];
				var nextv = (vv+1) % 3;
				var nextVert = verts[thisTri[nextv]];
				var displacement = [raystart[0]-thisVert[0],raystart[1]-thisVert[1],raystart[2]-thisVert[2]];
				var edgevector = [nextVert[0]-thisVert[0],nextVert[1]-thisVert[1],nextVert[2]-thisVert[2]];
				var crossProd = [ displacement[1]*edgevector[2] - displacement[2]*edgevector[1],
									displacement[2]*edgevector[0] - displacement[0]*edgevector[2],
									displacement[0]*edgevector[1] - displacement[1]*edgevector[0]];
				
				var dotProd = crossProd[0]*rayVec[0] + crossProd[1]*rayVec[1] + crossProd[2]*rayVec[2];
				
				//if (dotProd < 0){isInside=false;}
				signsSum+=dotProd/Math.abs(dotProd);
			}
			//isInsideTri[tt] = isInside;
			//if (isInside){
			if (Math.abs(signsSum)>2.5){	//should need 3, but unsure how numerical error works
				collidedFaces++;
				
				/*
				
				//draw simple collision point. 
				//calculate surface normal. 
				
				//plane defined by point-on-surface.normal = constant.
				//, can be expressed as Ax+By+Cy=D
				
				// ray:
				//raystart + t * rayDirection
				
				//solve for intersection:
				// norm.raystart + t*norm.rayDirection = D = norm.point-on-surface
				// =>
				// t = (norm.point-on-surface - norm.raystart)/norm.rayDirection
				// t = norm.(point-on-surface - raystart) / norm.rayDirection
				
				//collision point = 
				// raystart + rayDirection * (norm.(point-on-surface - raystart) / norm.rayDirection)
				
				//this is undefined if norm.rayDirection = 0 (ray is perpendicular to surface normal)
				// note that no need to normalise surface normal
				*/
				
				var vert0 = verts[thisTri[0]];
				var vert1 = verts[thisTri[1]]
				var vert2 = verts[thisTri[2]]
				
				/*
				var displacement = [raystart[0]-vert0[0],raystart[1]-vert0[1],raystart[2]-vert0[2]];
				
				*/
				//calc normal
				var dir1 = [ vert1[0]-vert0[0], vert1[1]-vert0[1], vert1[2]-vert0[2] ];
				var dir2 = [ vert2[0]-vert0[0], vert2[1]-vert0[1], vert2[2]-vert0[2] ];
				var faceNormal = [ dir1[1]*dir2[2] - dir1[2]*dir2[1],
									dir1[2]*dir2[0] - dir1[0]*dir2[2],
									dir1[0]*dir2[1] - dir1[1]*dir2[0]];
				
				//normalise face normal???
				var faceNormLenSq = faceNormal[0]*faceNormal[0] + faceNormal[1]*faceNormal[1] + faceNormal[2]*faceNormal[2]; 
				var faceNormInvLen = 1/Math.sqrt(faceNormLenSq);
				faceNormal[0]*=faceNormInvLen;
				faceNormal[1]*=faceNormInvLen;
				faceNormal[2]*=faceNormInvLen;
				
				/*
				var normDotRayDir = faceNormal[0]*rayVec[0]+ faceNormal[1]*rayVec[1]+ faceNormal[2]*rayVec[2];
				var normDotDisp = faceNormal[0]*displacement[0]+ faceNormal[1]*displacement[1]+ faceNormal[2]*displacement[2];
				var tval = -normDotDisp/normDotRayDir;	//- sign because point-on-surface calculated -ve	TODO check for 0
				var collisionPoint = [ raystart[0] + rayVec[0]*tval,
										raystart[1] + rayVec[1]*tval,
										raystart[2] + rayVec[2]*tval];
						//TODO WHAT IS THIS POINT?!!
				
				
				//gl.uniform4fv(shaderProgramColored.uniforms.uColor, [0.0, 1.0, 0.0, 1.0]);
//				drawSphere(collisionPoint);
		
				//simple ray/plane intersection maybe insufficient, because if ray is in plane, this point is undefined
				//do "closest approach" for space squashed about face normal - basically collision with the closest possible oblate ellipsoid centred on face. will result in bit odd glancing collision, but handle ray in plane case acceptably, without requiring special case.
				
				
				// displacement = raystart - faceCentre
				// stretchedDisplacement = displacement + displacement.normal * const * normal
				// stretchedRayvec = rayvec + rayvec.normal * const *normal
				
				//formulate in t (fraction of rayvec), move this fraction in unstretched space
				
				// t = stretchedDisplacement.stretchedRayvec / stretchedRayvec.stretchedRayvec (?)
				
				var normSq = faceNormal[0]*faceNormal[0] + faceNormal[1]*faceNormal[1] + faceNormal[2]*faceNormal[2];
				
				var const1 = 1000/normSq;	//a big number
				
				var dispDotNormTimesConst1 = const1 * ( displacement[0]*faceNormal[0]+ displacement[1]*faceNormal[1]+ displacement[2]*faceNormal[2] );
				
				var stretchedDisplacement = [ displacement[0] + dispDotNormTimesConst1*faceNormal[0],
												displacement[1] + dispDotNormTimesConst1*faceNormal[1],
												displacement[2] + dispDotNormTimesConst1*faceNormal[2]];
				
				var rayvecDotNormTimesConst1 = const1 * ( rayVec[0]*faceNormal[0]+ rayVec[1]*faceNormal[1]+ rayVec[2]*faceNormal[2] );
				
				var stretchedRayvec = [ rayVec[0] + rayvecDotNormTimesConst1*faceNormal[0],
										rayVec[1] + rayvecDotNormTimesConst1*faceNormal[1],
										rayVec[2] + rayvecDotNormTimesConst1*faceNormal[2]];
				
				var tt = ( stretchedDisplacement[0]*stretchedRayvec[0] + stretchedDisplacement[1]*stretchedRayvec[1] + stretchedDisplacement[2]*stretchedRayvec[2])/ ( stretchedRayvec[0]*stretchedRayvec[0] + stretchedRayvec[1]*stretchedRayvec[1] + stretchedRayvec[2]*stretchedRayvec[2]);
				 
				var colpoint = [ raystart[0] - tt*rayVec[0],
								raystart[1] - tt*rayVec[1],
								raystart[2] - tt*rayVec[2]];
						*/
				//collisionPoints.push(colpoint);
				//collisionPoints.push(true);	//just count a point. TODO check whether within line limits
				
				//to do this, for this face, calculate the normal, dot with difference between point on face and line start/end points. signs should be different (ie multiple should be negative)
				
				var dotOne = faceNormal[0]*(raystart[0]-vert0[0]) + faceNormal[1]*(raystart[1]-vert0[1]) + faceNormal[2]*(raystart[2]-vert0[2]);	//TODO simplifying with known start/end (eg fix to 0,0,0)
				var dotTwo = faceNormal[0]*(rayend[0]-vert0[0]) + faceNormal[1]*(rayend[1]-vert0[1]) + faceNormal[2]*(rayend[2]-vert0[2]);
					//TODO simplify by calculating facenormal dot vert0 separately
				if (dotOne*dotTwo > 0){
					collisionPoints.push(true);
				}
			}
			
		};
		return collisionPoints;	//TODO don't bother with points - just tell if collided		
	}
	
	return function(ii,jj,kk){
			//TODO should above stuff go inside an inner IIFE so unneeded stuff falls out of scope?
		var numcollisionpoints = checkForCollisions([(ii-32)/32,(jj-16)/32,(kk-32)/32],[0.02,0.5,0.01]).length;	//tiny extra 0.01,0.02 because suspect that rays through edges causing bad face count evenness
		return numcollisionpoints%2 == 0;
	}
})();


var seedValue;

function init(){
	
	stats = new Stats();
	stats.showPanel( 0 ); // 0: fps, 1: ms, 2: mb, 3+: custom
	document.body.appendChild( stats.dom );
	
	var gui = new dat.GUI();
	gui.add(guiParams, "box");
	gui.add(guiParams, "stupid");
	gui.add(guiParams, "sparse");
	gui.add(guiParams, "sparseWithNormals");
	gui.add(guiParams, "constrainToBox");
	
	canvas=document.getElementById("glcanvas");
	gl=glcanvas.getContext("webgl");
	
	
	
	canvas.addEventListener("mousedown", function(evt){
		mouseInfo.x = evt.offsetX;
		mouseInfo.y = evt.offsetY;
		mouseInfo.dragging = evt.buttons & 1;
		mouseInfo.lastPointingDir = getPointingDirectionFromScreenCoordinate({x:mouseInfo.x, y: mouseInfo.y});
		mouseInfo.buttons = evt.buttons;
		evt.preventDefault();
	});
	canvas.addEventListener("mouseup", function(evt){
		mouseInfo.dragging = evt.buttons & 1;
		mouseInfo.buttons = evt.buttons;
	});
	canvas.addEventListener("mouseout", function(evt){
		mouseInfo.dragging = false;
		mouseInfo.buttons = 0;
	});
	canvas.addEventListener("mousemove", function(evt){
		mouseInfo.currentPointingDir = getPointingDirectionFromScreenCoordinate({x:evt.offsetX, y: evt.offsetY});
		if (mouseInfo.dragging){
			var pointingDir = mouseInfo.currentPointingDir;
			//console.log("pointingDir = " + pointingDir);
			
			//get the direction of current and previous mouse position.
			//do a cross product to work out the angle rotated
			//and rotate the player by this amount
			
			var crossProd = crossProductHomgenous(pointingDir, mouseInfo.lastPointingDir);
			mouseInfo.lastPointingDir = pointingDir;
			
			//rotate player 
			//guess have signs here because of unplanned handedness of screen, 3d co-ord systems
			var rotateAmt = [-crossProd.x / crossProd.w, crossProd.y / crossProd.w, crossProd.z / crossProd.w];
			rotatePlayer(rotateAmt);
			
		}
		/*
		if (pointerLocked){
			rotatePlayer([ -0.001* evt.movementY, -0.001* evt.movementX, 0]);	//TODO screen resolution dependent sensitivity.
		}
		*/
	});
	canvas.addEventListener("click", function(evt){
		//movePlayer([0,0,0.05]);	//basic click to move forward. TODO keyboard controls
	});
	
	
	crosssectionscanvas=document.getElementById("mycanvas");
	ctx=crosssectionscanvas.getContext("2d");
	
	
	crosssectionscanvas.width = canvassize;
	crosssectionscanvas.height = canvassize;
	
	//voxdata = [];
	//generate arrays 
	for (var ii=0;ii<blocksize;ii++){
		var slicedata=[];
		voxdata.push(slicedata);
		for (var jj=0;jj<blocksize;jj++){
			slicedata.push([]);	//could push new Array(blocksize). use empty array if fill with sparse data. unimportant for small blocksize 
		}
	}
	
	//fill with data
	function makeVoxdataForFunc(thefunction){
		for (var ii=0;ii<blocksize;ii++){
			var slicedata = voxdata[ii];
			for (var jj=0;jj<blocksize;jj++){
				var stripdata = slicedata[jj];
				for (var kk=0;kk<blocksize;kk++){
					stripdata[kk] = (thefunction(ii,jj,kk) > 0 );	//if know starting from empty array, can use push. if want sparse array, should use conditional. 
				}											//note some functions may benefit from currying eg functionForXY= thefunction(ii,jj), 
			}
		}
	}
	
	//var voxFunction = sinesfunction;
	//var voxFunction = landscapeFunction;
	//var voxFunction = bigBallFunction;
	//var voxFunction = bigCylinderFunction;
	//var voxFunction = curveCornerFunction;
	//var voxFunction = bilinearFilterBinaryFunctionGen(sinesfunction);
	//var voxFunction = bilinearFilterBinaryFunctionGen(landscapeFunction);
	//var voxFunction = bilinearFilterBinaryFunctionGen(bigBallFunction);
	//var voxFunction = bilinearFilterBinaryFunctionGen(bigCylinderFunction);
	//var voxFunction = bilinearFilterBinaryFunctionGen(curveCornerFunction);
	var voxFunction = bilinearFilterBinaryFunctionGen(perlinPlanetFunction);
	//var voxFunction = perlinPlanetFunction;
	
	seedValue= Math.random();
	console.log("seed: " + seedValue);
	var genStartTime = Date.now();
	//noise.seed(seedValue);var voxFunction = perlinfunction;
	
	//var voxFunction = fromPolyModelFunction;
	var voxFunction = fromPolyModelFunctionFast;
	
	makeVoxdataForFunc(voxFunction);	
	console.log("Time taken to generate: " + (Date.now()-genStartTime));
	
	function perlinfunction(ii,jj,kk){
		//return 10*noise.perlin3(ii/12,jj/12,kk/12);	//if divide by too small number, too many indices generated
		return 10*noise.perlin3(ii/24,jj/12,kk/12) - 0.75*kk +20;	//landscape with 3d perlin surface
	}
	function perlinPlanetFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		kkm = kk-blocksize/2;
		//return (800+200*noise.perlin3(ii/12,jj/12,kk/12)) - (iim*iim + jjm*jjm + kkm*kkm);	//lumpy
		//return (500+500*noise.perlin3(ii/6,jj/6,kk/6)) - (iim*iim + jjm*jjm + kkm*kkm);
		//return (500+500*noise.perlin3(ii/12,jj/12,kk/12)) - (iim*iim + jjm*jjm + kkm*kkm);	//blobby
		return Math.min(500, 500+1000*noise.perlin3(ii/6,jj/6,kk/6)) - (iim*iim + jjm*jjm + kkm*kkm);	//not good because discontinuities affect gradient
		//return (500+500*noise.perlin3(ii/12,jj/12,kk/12)) - 0.001*Math.pow((iim*iim + jjm*jjm + kkm*kkm),2);	//not good because discontinuities affect gradient
	}
	
	function sinesfunction(ii,jj,kk){
		var sinscale=3;
		//return Math.sin(ii/sinscale)+Math.sin(jj/sinscale)+Math.sin(kk/sinscale);
		return Math.sin(ii/sinscale)+Math.sin(jj/sinscale)+Math.sin(kk/sinscale) - (ii/20 - 0.5);
	}
	function bigBallFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		kkm = kk-blocksize/2;
		return iim*iim + jjm*jjm + kkm*kkm - 1000;
	}
	function bigCylinderFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		kkm = kk-blocksize/2;
		//return Math.max(iim*iim + jjm*jjm, kkm*kkm) - 1000;	//inverted
		//return 950- Math.max(iim*iim + jjm*jjm, kkm*kkm);
		return 950- Math.max( 0.55*(iim*iim + jjm*jjm + kkm*kkm), Math.max(iim*iim + jjm*jjm, kkm*kkm)); //bevel
	}
	function curveCornerFunction(ii,jj,kk){
		return 10000-ii*jj*kk;
	}
	function landscapeFunction(ii,jj,kk){	
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		kkm = kk-blocksize/2 +15;	//add to move hole deeper
		var ballPot = Math.sqrt(iim*iim + jjm*jjm + kkm*kkm) - 20;
		
		var sinscale=3;
		var eggboxPot = Math.sin(ii/sinscale)+Math.sin(jj/sinscale);
		
		//return Math.min( eggboxPot - (kk/2 - 15) , ballPot/2);	//sharp field causes jagged edge when using standard "downhill" vertex adjustment (which assumes constant gradient)
		
		return ((eggboxPot - (kk/2 - 15))+3) * ((ballPot/2)+3) - 9;	//smoother field
		
	}
	
	function bilinearFilterBinaryFunctionGen(smoothFunction){	//generate a function that returns 1,-1 for occupied/unoccupied grid points, bilinear smoothed value between
		//TODO pregenerate all data (64x64x64) inside function - saves doing 8 function evaluations every time run the returned function
		return function(ii,jj,kk){
			iif=Math.floor(ii);
			jjf=Math.floor(jj);
			kkf=Math.floor(kk);
			//return myVoxData[iif][jjf][kkf];	//todo bilinear filter
			
			//remainder
			var iir = ii-iif;
			var jjr = jj-jjf;
			var kkr = kk-kkf;
			
			//bilinear filter version (bit slower. could code dedicated gradient extraction more efficiently)
			if (iif<63 && jjf<63 & kkf<63){
				var sum=0;
				/*
				sum+= myVoxData[iif][jjf][kkf]*(1-iir)*(1-jjr)*(1-kkr);
				sum+= myVoxData[iif][jjf][kkf+1]*(1-iir)*(1-jjr)*kkr;
				sum+= myVoxData[iif][jjf+1][kkf]*(1-iir)*jjr*(1-kkr);
				sum+= myVoxData[iif][jjf+1][kkf+1]*(1-iir)*jjr*kkr;
				sum+= myVoxData[iif+1][jjf][kkf]*iir*(1-jjr)*(1-kkr);
				sum+= myVoxData[iif+1][jjf][kkf+1]*iir*(1-jjr)*kkr;
				sum+= myVoxData[iif+1][jjf+1][kkf]*iir*jjr*(1-kkr);
				sum+= myVoxData[iif+1][jjf+1][kkf+1]*iir*jjr*kkr;
				*/
				sum+= (smoothFunction(iif,jjf,kkf) <0 ? -1:1)*(1-iir)*(1-jjr)*(1-kkr);
				sum+= (smoothFunction(iif,jjf,kkf+1) <0 ? -1:1) *(1-iir)*(1-jjr)*kkr ;
				sum+= (smoothFunction(iif,jjf+1,kkf) <0 ? -1:1) *(1-iir)*jjr*(1-kkr);
				sum+= (smoothFunction(iif,jjf+1,kkf+1) <0 ? -1:1)*(1-iir)*jjr*kkr;
				sum+= (smoothFunction(iif+1,jjf,kkf) <0 ? -1:1)*iir*(1-jjr)*(1-kkr) ;
				sum+= (smoothFunction(iif+1,jjf,kkf+1) <0 ? -1:1)*iir*(1-jjr)*kkr;
				sum+= (smoothFunction(iif+1,jjf+1,kkf) <0 ? -1:1)*iir*jjr*(1-kkr);
				sum+= (smoothFunction(iif+1,jjf+1,kkf+1) <0? -1:1)*iir*jjr*kkr;
				return sum;
			}else{
				//strictly there are more edge cases but don't care for now
				//maybe this never gets hit... 
				//return myVoxData[63][63][63];
				return smoothFunction(63,63,63);
				console.log("EDGE CASE HIT");
			}
			
		}
	}
	
	console.log(voxdata);
	
	//draw to canvas to test algo.
	//start with a single slice
	var imgData = ctx.createImageData(blocksize, blocksize);
	var imgDataData = imgData.data;
	for (var ii=0;ii<blocksize;ii++){
		var slicedata = voxdata[ii];
		var idx=0;
		for (var jj=0;jj<blocksize;jj++){
			var stripdata = slicedata[jj];
			for (var kk=0;kk<blocksize;kk++){
				var color = stripdata[kk] ? 255 : 0;
				imgDataData[idx]=color;
				imgDataData[idx+1]=color;
				imgDataData[idx+2]=color;
				imgDataData[idx+3]=255;
				idx+=4;
			}
		}
		ctx.putImageData(imgData, blocksize*(ii%8), blocksize*((ii - ii%8)/8));
	}
	
	
	//generate mesh data for voxels. TODO put this into another file
	//create (sparse) array containing the vertex id for each 3d point
	//either as 2nd step or as go along, create 2 tris between any pairs of occupied/unoccupied blocks
	
	//try just counting faces. ignore outermost
	//if having dedicated verts for faces in each direction, quads<=verts<=4*quads	. for roundish objects, suspect closer to latter
	var totalQuads=0;
	//TODO 3 sets of faces
	for (var ii=0;ii<blocksize-1;ii++){
		var slicedata = voxdata[ii];
		var nextSlicedata = voxdata[ii+1];
		for (var jj=0;jj<blocksize;jj++){
			var stripData=slicedata[jj];
			var nextStripData=nextSlicedata[jj];
			for (var kk=0;kk<blocksize;kk++){
				if (stripData[kk]!=nextStripData[kk]){
					totalQuads++;
				}
			}
		}
	}
	console.log(totalQuads);
	for (var ii=0;ii<blocksize;ii++){
		var slicedata = voxdata[ii];
		for (var jj=0;jj<blocksize-1;jj++){
			var stripData=slicedata[jj];
			var nextStripData=slicedata[jj+1];
			for (var kk=0;kk<blocksize;kk++){
				if (stripData[kk]!=nextStripData[kk]){
					totalQuads++;
				}
			}
		}
	}
	console.log(totalQuads);
	for (var ii=0;ii<blocksize;ii++){
		var slicedata = voxdata[ii];
		for (var jj=0;jj<blocksize;jj++){
			var stripData=slicedata[jj];
			for (var kk=0;kk<blocksize-1;kk++){
				if (stripData[kk]!=stripData[kk+1]){
					totalQuads++;
				}
			}
		}
	}
	console.log(totalQuads);
	//~49k for array of spheres.
	//~18k for big ball
	
	//todo count unique vertices (in all 6 directions)
	
	//todo not indexed version (easier, no problem with 65536 index limit)
	
	voxData["stupid"] = (function(){
		//stupid implementation. a vertex for every grid point, regardless of whether occupied.
		//can only do part of 64x64x64 this way
		var vertices = [];
		var indices = [];
		
		for (var ii=0;ii<=32;ii++){
			for (var jj=0;jj<=32;jj++){
				for (var kk=0;kk<=32;kk++){
					vertices.push(ii/32, jj/32, kk/32);
				}
			}
		}
		
		var oneVertIdx;
		for (var ii=0;ii<32;ii++){
			for (var jj=0;jj<32;jj++){
				for (var kk=0;kk<31;kk++){
					var difference=voxdata[ii][jj][kk+1]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = ii*33*33 + jj*33 + kk + 1 ;
						if ( difference>0 ){
							indices.push( oneVertIdx , oneVertIdx+33 , oneVertIdx+33+33*33 );	//bottom faces
							indices.push( oneVertIdx , oneVertIdx+33+33*33, oneVertIdx+33*33);
						}else{
							indices.push( oneVertIdx, oneVertIdx+33+33*33, oneVertIdx+33 );		//top faces
							indices.push( oneVertIdx, oneVertIdx+33*33, oneVertIdx+33+33*33);
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<32;ii++){
			for (var jj=0;jj<31;jj++){
				for (var kk=0;kk<32;kk++){
					var difference=voxdata[ii][jj+1][kk]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = ii*33*33 + jj*33 + kk + 33 ;
						if ( difference<0 ){
							indices.push( oneVertIdx, oneVertIdx+1, oneVertIdx+1+33*33 );
							indices.push( oneVertIdx, oneVertIdx+1+33*33, oneVertIdx+33*33);
						}else{
							indices.push( oneVertIdx, oneVertIdx+1+33*33, oneVertIdx+1 );
							indices.push( oneVertIdx, oneVertIdx+33*33, oneVertIdx+1+33*33);
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<31;ii++){
			for (var jj=0;jj<32;jj++){
				for (var kk=0;kk<32;kk++){
					var difference=voxdata[ii+1][jj][kk]-voxdata[ii][jj][kk];
					if ( difference!=0 ){
						oneVertIdx = ii*33*33 + jj*33 + kk + 33*33 ;					
						if ( difference>0 ){
							indices.push( oneVertIdx, oneVertIdx+1, oneVertIdx+1+33);
							indices.push( oneVertIdx, oneVertIdx+1+33, oneVertIdx+33);
						}else{
							indices.push( oneVertIdx, oneVertIdx+1+33, oneVertIdx+1);
							indices.push( oneVertIdx, oneVertIdx+33, oneVertIdx+1+33);
						}
					}
				}
			}
		}
		
		return {
			vertices:vertices,
			indices:indices
		};
	})();	
	
	var topFaceCount=0;
	voxData["sparse"] = (function(){
		//stupid implementation. a vertex for every grid point, regardless of whether occupied.
		//can only do part of 64x64x64 this way
		var vertices = [];
		var indices = [];
		var normals = [];
		var colors = [];
		
		//sparse version - no unused vertices.
		var indexForGridPoint = [];
		var currentPoint = 0;
		var numNeighbours;
		var occNeighbours;
		for (var ii=0;ii<=64;ii++){
			for (var jj=0;jj<=64;jj++){
				for (var kk=0;kk<=64;kk++){
					numNeighbours=1;
					occNeighbours=0;
					if (ii>0 && ii<64){numNeighbours*=2;};
					if (jj>0 && jj<64){numNeighbours*=2;};
					if (kk>0 && kk<64){numNeighbours*=2;};
					
					if (ii<64 && jj<64 && kk<64 && voxdata[ii][jj][kk]){
						occNeighbours++;
					}
					if (ii<64 && jj<64 && kk>0 && voxdata[ii][jj][kk-1]){
						occNeighbours++;
					}
					if (ii<64 && jj>0 && kk<64 && voxdata[ii][jj-1][kk]){
						occNeighbours++;
					}
					if (ii<64 && jj>0 && kk>0 && voxdata[ii][jj-1][kk-1]){
						occNeighbours++;
					}
					if (ii>0 && jj<64 && kk<64 && voxdata[ii-1][jj][kk]){
						occNeighbours++;
					}
					if (ii>0 && jj<64 && kk>0 && voxdata[ii-1][jj][kk-1]){
						occNeighbours++;
					}
					if (ii>0 && jj>0 && kk<64 && voxdata[ii-1][jj-1][kk]){
						occNeighbours++;
					}
					if (ii>0 && jj>0 && kk>0 && voxdata[ii-1][jj-1][kk-1]){
						occNeighbours++;
					}
					if (occNeighbours%numNeighbours){	// ( !=0 )
						indexForGridPoint[getNumberOfGridPoint(ii,jj,kk)] = currentPoint;
						addVertData(ii,jj,kk);
						currentPoint++;
					}	
				}
			}
		}
		console.log(currentPoint);
		console.log(indexForGridPoint);
		
		function addVertData(ii,jj,kk){
			ii-=0.5;	//???
			jj-=0.5;
			kk-=0.5;
			
			//vertices.push(ii/32, jj/32, kk/32);
			//normals.push(0,0,0);
				//^^ little faster for doing O(3) slow teapot thing

			//smooth vertices (TODO make 2 vert buffers and UI to switch between) 
			//just look at gradient between surrounding points, move downhill. or take analytic gradient from something function used to generate vox data
			//to make this work generally without needing to calculate analytic derivatives, just use numerical differences.
			var delta = 0.01;
			var fudgeFactor = 0.7;	//less than 1 to avoid overshoot
			var centralPoint,gradX,gradY,gradZ,totalGradSq,sharedPart;
			
			for (var iter=0;iter<3;iter++){
				centralPoint = voxFunction(ii,jj,kk);
							
				gradX = (voxFunction(ii+delta,jj,kk)- centralPoint)/delta;
				gradY = (voxFunction(ii,jj+delta,kk)- centralPoint)/delta;
				gradZ = (voxFunction(ii,jj,kk+delta)- centralPoint)/delta;
				
				totalGradSq = gradX*gradX + gradY*gradY + gradZ*gradZ;
				
				//have some sort of hill normal. should move downhill
				//move by something like (discrepancy / totalGradent)*(gradientVector/totalGradient)
				//to avoid /0 error add something to totalGradient
			
				sharedPart = centralPoint / ( totalGradSq + 0.001);
				ii = ii-sharedPart*gradX*fudgeFactor;
				jj = jj-sharedPart*gradY*fudgeFactor;
				kk = kk-sharedPart*gradZ*fudgeFactor;
			}
			
			vertices.push(ii/32, jj/32, kk/32);
			var invLength = Math.sqrt(1/( totalGradSq + 0.001));
			normals.push(invLength*gradX, invLength*gradY, invLength*gradZ);
			
			//colors.push(Math.random(),Math.random(),Math.random());
			var colorScale = 6;	//scale of noise (smaller = finer)
			//var grayColor = 0.5+0.5*noise.perlin3(ii/colorScale,jj/colorScale,kk/colorScale);	// mapt -1 to 1 -> 0 to 1
			var grayColor = 1.5+0.5*noise.perlin3(ii/colorScale,jj/colorScale,kk/colorScale);
			colors.push(grayColor * 0.5,
						grayColor * 0.05,
						grayColor * 0.01
			);
			//normals.push(0,0,0); 
			
		}
		function getNumberOfGridPoint(ii,jj,kk){
			return ii*65*65 + jj*65 + kk;
		}
		function pushIndexForNumber(nii,njj,nkk){
			indices.push(indexForGridPoint[nii],indexForGridPoint[njj],indexForGridPoint[nkk]);
		}
		
		var oneVertIdx;
		for (var ii=0;ii<64;ii++){
			for (var jj=0;jj<64;jj++){
				for (var kk=0;kk<63;kk++){
					var difference=voxdata[ii][jj][kk+1]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = getNumberOfGridPoint(ii,jj,kk+1);
						if ( difference>0 ){
							pushIndexForNumber( oneVertIdx , oneVertIdx+65 , oneVertIdx+65+65*65 );	//bottom faces
							pushIndexForNumber( oneVertIdx , oneVertIdx+65+65*65, oneVertIdx+65*65);
						}else{
							pushIndexForNumber( oneVertIdx, oneVertIdx+65+65*65, oneVertIdx+65 );		//top faces
							pushIndexForNumber( oneVertIdx, oneVertIdx+65*65, oneVertIdx+65+65*65);
							topFaceCount++;
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<64;ii++){
			for (var jj=0;jj<63;jj++){
				for (var kk=0;kk<64;kk++){
					var difference=voxdata[ii][jj+1][kk]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = getNumberOfGridPoint(ii,jj+1,kk);
						if ( difference<0 ){
							pushIndexForNumber( oneVertIdx, oneVertIdx+1, oneVertIdx+1+65*65 );
							pushIndexForNumber( oneVertIdx, oneVertIdx+1+65*65, oneVertIdx+65*65);
						}else{
							pushIndexForNumber( oneVertIdx, oneVertIdx+1+65*65, oneVertIdx+1 );
							pushIndexForNumber( oneVertIdx, oneVertIdx+65*65, oneVertIdx+1+65*65);
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<63;ii++){
			for (var jj=0;jj<64;jj++){
				for (var kk=0;kk<64;kk++){
					var difference=voxdata[ii+1][jj][kk]-voxdata[ii][jj][kk];
					if ( difference!=0 ){
						oneVertIdx = getNumberOfGridPoint(ii+1,jj,kk);					
						if ( difference>0 ){
							pushIndexForNumber( oneVertIdx, oneVertIdx+1, oneVertIdx+1+65);
							pushIndexForNumber( oneVertIdx, oneVertIdx+1+65, oneVertIdx+65);
						}else{
							pushIndexForNumber( oneVertIdx, oneVertIdx+1+65, oneVertIdx+1);
							pushIndexForNumber( oneVertIdx, oneVertIdx+65, oneVertIdx+1+65);
						}
					}
				}
			}
		}
		
		return {
			vertices:vertices,
			normals:normals,
			colors:colors,
			indices:indices
		};
	})();
	
	console.log("topFaceCount: " + topFaceCount);
	
	
	initGL();
	
	gl.clearColor.apply(gl,[1,1,0,1]);
	
	initShaders();
	initBuffers();
	
	//hack - draw with most complex shader prog. this binds buffers for all used attribute indices
	//possibly, better to handle which are enabled when switching shader using enableVertexAttribArray / disableVertexAttribArray
	//drawing with unneeded attribArrays bound still works
	gl.useProgram(shaderPrograms.withNormalsAndColor);
	drawObjectFromBuffers(voxBuffers["sparseWithNormals"], shaderPrograms.withNormalsAndColor);
	
	gl.enable(gl.DEPTH_TEST);
	
	requestAnimationFrame(itMechanicsAndDrawScene);
}

//copied from 3sphere project. much of this unused since objs here don't have normals, textures
function drawObjectFromBuffers(bufferObj, shaderProg, usesCubeMap){
	prepBuffersForDrawing(bufferObj, shaderProg, usesCubeMap);
	drawObjectFromPreppedBuffers(bufferObj, shaderProg);
}
function prepBuffersForDrawing(bufferObj, shaderProg, usesCubeMap){
									
	gl.enable(gl.CULL_FACE);
	gl.bindBuffer(gl.ARRAY_BUFFER, bufferObj.vertexPositionBuffer);
    gl.vertexAttribPointer(shaderProg.attributes.aVertexPosition, bufferObj.vertexPositionBuffer.itemSize, gl.FLOAT, false, 0, 0);
	
	if (bufferObj.vertexNormalBuffer && shaderProg.attributes.aVertexNormal){
		gl.bindBuffer(gl.ARRAY_BUFFER, bufferObj.vertexNormalBuffer);
		gl.vertexAttribPointer(shaderProg.attributes.aVertexNormal, bufferObj.vertexNormalBuffer.itemSize, gl.FLOAT, false, 0, 0);
	}
	if (bufferObj.vertexColorBuffer && shaderProg.attributes.aVertexColor){
		gl.bindBuffer(gl.ARRAY_BUFFER, bufferObj.vertexColorBuffer);
		gl.vertexAttribPointer(shaderProg.attributes.aVertexColor, bufferObj.vertexColorBuffer.itemSize, gl.FLOAT, false, 0, 0);
	}
	
	gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, bufferObj.vertexIndexBuffer);
	/*
	if (bufferObj.vertexTextureCoordBuffer){
		gl.bindBuffer(gl.ARRAY_BUFFER, bufferObj.vertexTextureCoordBuffer);
		gl.vertexAttribPointer(shaderProg.attributes.aTextureCoord, bufferObj.vertexTextureCoordBuffer.itemSize, gl.FLOAT, false, 0, 0);
		gl.activeTexture(gl.TEXTURE0);
		gl.bindTexture(gl.TEXTURE_2D, texture);
		gl.uniform1i(shaderProg.uniforms.uSampler, 0);
	}
	
	if (usesCubeMap){
		gl.activeTexture(gl.TEXTURE0);
		gl.bindTexture(gl.TEXTURE_CUBE_MAP, cubemapTexture);
		gl.uniform1i(shaderProg.uniforms.uSampler, 0);
	}
	*/
	//gl.uniformMatrix4fv(shaderProg.uniforms.uPMatrix, false, pMatrix);
}


function drawObjectFromPreppedBuffers(bufferObj, shaderProg){
	gl.uniformMatrix4fv(shaderProg.uniforms.uPMatrix, false, pMatrix);	//TODO not every frame.
	gl.uniformMatrix4fv(shaderProg.uniforms.uMVMatrix, false, mvMatrix);
	gl.drawElements(gl.TRIANGLES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
}

var currentTime=0;

//borrowed from rayMarcherTest
var iterateMechanics = (function generateIterateMechanics(){
	var newTime = Date.now();
	var camSpeed = [0,0,0];
	
	window.addEventListener("keydown",function(e){
		var keyCode=e.keyCode;
		console.log(keyCode);
		//WASD = 87,65,83,68
		var camStep = 0.001;
		
		var camMove = [0,0,0];
		if (keyCode == 87){camMove[2]+=camStep;}
		if (keyCode == 83){camMove[2]-=camStep;}
		if (keyCode == 65){camMove[0]+=camStep;}
		if (keyCode == 68){camMove[0]-=camStep;}
		
		//add to camSpeed instead
		/*
		camSpeed[0] += camMove[0]*mvMatrix[0] + camMove[1]*mvMatrix[4] + + camMove[2]*mvMatrix[8];
		camSpeed[1] += camMove[0]*mvMatrix[1] + camMove[1]*mvMatrix[5] + + camMove[2]*mvMatrix[9];
		camSpeed[2] += camMove[0]*mvMatrix[2] + camMove[1]*mvMatrix[6] + + camMove[2]*mvMatrix[10];
		*/
		//transposed vs code from rayMarcherTest
		camSpeed[0] += camMove[0]*mvMatrix[0] + camMove[1]*mvMatrix[1] + + camMove[2]*mvMatrix[2];
		camSpeed[1] += camMove[0]*mvMatrix[4] + camMove[1]*mvMatrix[5] + + camMove[2]*mvMatrix[6];
		camSpeed[2] += camMove[0]*mvMatrix[8] + camMove[1]*mvMatrix[9] + + camMove[2]*mvMatrix[10];
		
		var rollSpd = 0.05;
		var rollAmt = 0;
		if (keyCode == 69){rollAmt-=rollSpd;}
		if (keyCode == 81){rollAmt+=rollSpd;}
		rotatePlayer([0,0,-rollAmt]);	// - because suspect matrix inverted compared to rayMarcherTest...
	});
	
	return function(){
		var oldTime = newTime;
		newTime = Date.now();
		var timeElapsed = Math.min(newTime - oldTime, 1000);	//1s max
		//exponential decay of speed towards desired speed. TODO proper movement integration, but doesn't really matter
		var decayFactor = Math.pow(0.995,timeElapsed);
		
		camSpeed = camSpeed.map(function(val,ii){return val*decayFactor;});	//TODO keystate = target value... + (1-decayFactor)*camMove[ii];});
		//camSpeed = camSpeed.map(function(val,ii){return camMove[ii];});
		playerPosition = playerPosition.map(function(val,ii){return val+timeElapsed*camSpeed[ii];});
		//playerPosition = playerPosition.map(function(val,ii){return val;});
		
		if (guiParams.constrainToBox){
			for (var cc=0;cc<3;cc++){
				playerPosition[cc] = Math.max(-2.0,Math.min(0.0,playerPosition[cc]));	//constraint values wierd. guess messed up matrix!
			}
		}
		
		setPlayerTranslation(playerPosition);
	}
})();
	
function itMechanicsAndDrawScene(){
	requestAnimationFrame(itMechanicsAndDrawScene);
	iterateMechanics();
	drawScene();
	
	document.getElementById("debugtext").innerHTML = checkCollision(playerPosition);
}

function drawScene(drawTime){
	//requestAnimationFrame(drawScene);
	
	stats.end();
	stats.begin();
	
	resizecanvas(1);	//TODO should this really happen every frame? perf impact?
	gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);	//TODO should this be every frame?

	mat4.perspective(mainCamFov, gl.viewportWidth/gl.viewportHeight, 0.01,20.0,pMatrix);	//apparently 0.9.5, last param is matrix rather than 1st!! todo use newer???
																	//also old one uses degs!
																	
	mat4.set(playerMatrix, mvMatrix);																
																	
	gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
	
	var activeShaderProgram = shaderPrograms.simple;
	gl.useProgram(activeShaderProgram);
	
	//todo if drawing many, prep buffers once.
	if (guiParams.box){
		drawObjectFromBuffers(basicCubeBuffers, activeShaderProgram);
	}
	if (guiParams.stupid){
		drawObjectFromBuffers(voxBuffers["stupid"], activeShaderProgram);
	}
	if (guiParams.sparse){
		drawObjectFromBuffers(voxBuffers["sparse"], activeShaderProgram);
	}
	
	if (guiParams.sparseWithNormals){
		//activeShaderProgram = shaderPrograms.withNormals;
		activeShaderProgram = shaderPrograms.withNormalsAndColor;
		gl.useProgram(activeShaderProgram);
		drawObjectFromBuffers(voxBuffers["sparseWithNormals"], activeShaderProgram);
	}
	
	//mat4.rotateZ(mvMatrix,0.0003*(drawTime-currentTime)); 
	currentTime=drawTime;
}


function movePlayer(vec){	//[left,up,forward]
	playerPosition[0] += vec[0]*playerMatrix[0] + vec[1]*playerMatrix[1] + vec[2]*playerMatrix[2];
	playerPosition[1] += vec[0]*playerMatrix[4] + vec[1]*playerMatrix[5] + vec[2]*playerMatrix[6];
	playerPosition[2] += vec[0]*playerMatrix[8] + vec[1]*playerMatrix[9] + vec[2]*playerMatrix[10];
	setPlayerTranslation(playerPosition);
	//constrainPlayerPositionToBox();
}

function rotatePlayer(vec){
	setPlayerTranslation([0,0,0]);
	var rotationMag = Math.sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
	var rotMat=mat4.identity();
	mat4.rotate(rotMat, rotationMag, [vec[0]/rotationMag, vec[1]/rotationMag, vec[2]/rotationMag]);	//TODO find or make method taking vector instead of separate unit axis/angle
	mat4.multiply(rotMat, playerMatrix, playerMatrix);
	setPlayerTranslation(playerPosition);
}

function setPlayerTranslation(posArray){
	setMatTranslation(playerMatrix, posArray);
}

function setMatTranslation(mat, posArray){
	//probably using mat4 wierdly so standard mat4.translate isn't working sensibly
	mat[12]=0;mat[13]=0;mat[14]=0;
	mat4.translate(mat, posArray);
}

function getPointingDirectionFromScreenCoordinate(coords){
	
	var maxyvert = 1.0;	
	var maxxvert = screenAspect;
	
	var xpos = maxxvert*(coords.x*2.0/gl.viewportWidth   -1.0 );
	var ypos = maxyvert*(coords.y*2.0/gl.viewportHeight   -1.0 );
	var radsq = xpos*xpos + ypos*ypos;
	var zpos = 1.0/Math.tan(mainCamFov*Math.PI/360); //TODO precalc

	//normalise - use sending back homogenous co-ords because maybe a tiny amount more efficient since cross producting anyway
	var mag= Math.sqrt(radsq + zpos*zpos);
	
	return {
		x: xpos,
		y: ypos,
		z: zpos,
		w: mag
	}
}

function crossProductHomgenous(dir1, dir2){
	var output ={};
	output.x = dir1.y * dir2.z - dir1.z * dir2.y; 
	output.y = dir1.z * dir2.x - dir1.x * dir2.z; 
	output.z = dir1.x * dir2.y - dir1.y * dir2.x;
	output.w = dir1.w * dir2.w;
	return output;
}

function checkCollision(position){
	//look up what grid square are in, return voxdata for that. (TODO use smoothed/polygon data)
	//TODO check offset of 0.5 etc...
	//TODO make position sensible! 
	
	var ii = Math.floor(-position[0]*32);
	var jj = Math.floor(-position[1]*32);
	var kk = Math.floor(-position[2]*32);
	if (ii<0 || jj<0 || kk<0 || ii>63 || jj >63 || kk>63){
		//return false;
		return "OUTSIDE BOX. ii=" + ii + ", jj=" + jj + ", kk=" +kk;
	}
	return voxdata[ii][jj][kk];
}