
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
	voxBuffers["sparseSmoothed"]={};
	voxBuffers["sparseWithNormals"]={};
	voxBuffers["sparseBasicAvgSmoothed"]={};
	voxBuffers["sparseDCSmoothed"]={};
	voxBuffers["sparseDCWithNormals"]={};
	loadBufferData(voxBuffers["stupid"], voxData["stupid"]);
		
	loadBufferData(voxBuffers["sparse"], { vertices:voxData["sparse"].vertices, indices:voxData["sparse"].indices, directionalIndices:voxData["sparse"].directionalIndices});	//copy all but normals!
	loadBufferData(voxBuffers["sparseSmoothed"], { vertices:voxData["sparse"].smoothVertices, indices:voxData["sparse"].indices});	//copy all but normals!
	loadBufferData(voxBuffers["sparseWithNormals"], { vertices:voxData["sparse"].smoothVertices, indices:voxData["sparse"].indices, normals:voxData["sparse"].normals, colors:voxData["sparse"].colors});
	loadBufferData(voxBuffers["sparseBasicAvgSmoothed"], { vertices:voxData["sparse"].basicAvgVertices, indices:voxData["sparse"].indices});	//copy all but normals!
	loadBufferData(voxBuffers["sparseDCSmoothed"], { vertices:voxData["sparse"].dcVertices, indices:voxData["sparse"].indices});	//copy all but normals!
	loadBufferData(voxBuffers["sparseDCWithNormals"], { vertices:voxData["sparse"].dcVertices, indices:voxData["sparse"].indices, normals:voxData["sparse"].dcNormals, colors:voxData["sparse"].dcColors});
	
	//note could share buffers for some of above - currently generate multiple buffers from the same data
	
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
		
		if (sourceData.directionalIndices){
			bufferObj.directionalIndices = sourceData.directionalIndices;
		}
	}
}


var shaderPrograms={};
function initShaders(){
	shaderPrograms.simple = loadShader( "shader-simple-vs", "shader-simple-fs",{
					attributes:["aVertexPosition"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
					
	shaderPrograms.simpleColor = loadShader( "shader-simple-vs", "shader-simple-color-fs",{
					attributes:["aVertexPosition"],
					uniforms:["uMVMatrix","uPMatrix","uColor"]
					});
					
	shaderPrograms.withNormals = loadShader( "shader-withnormals-vs", "shader-color-fs",{
					attributes:["aVertexPosition","aVertexNormal"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
	
	shaderPrograms.withNormalsAndColor = loadShader( "shader-withnormals-andcolor-vs", "shader-color-fs",{
					attributes:["aVertexPosition","aVertexNormal","aVertexColor"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
					
	shaderPrograms.texmap = loadShader( "shader-texmap-vs", "shader-texmap-fs",{
					attributes:["aVertexPosition","aVertexNormal","aVertexColor"],
					uniforms:["uMVMatrix","uPMatrix"]
					});
	
	shaderPrograms.nmap = loadShader( "shader-nmap-vs", "shader-nmap-fs",{
					attributes:["aVertexPosition","aVertexNormal","aVertexColor"],
					uniforms:["uMVMatrix","uPMatrix"]
					});	
}

var texture;
function initTexture(){
	//texture = makeTexture("img/4483-v7.jpg");
	//texture = makeTexture("img/cretish0958.png");
	texture = makeTexture("img/normal_mapping_normal_map.png");	//brick
}
function makeTexture(src) {	//to do OO
	var texture = gl.createTexture();
	texture.image = new Image();
	texture.image.onload = function(){
		bind2dTextureIfRequired(texture);
		gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
		
		gl.pixelStorei(gl.UNPACK_COLORSPACE_CONVERSION_WEBGL, gl.NONE);	//linear colorspace grad light texture (TODO handle other texture differently?)
		
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_LINEAR);
		//gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
		gl.generateMipmap(gl.TEXTURE_2D);
		bind2dTextureIfRequired(null);
	};	
	texture.image.src = src;
	return texture;
}
var bind2dTextureIfRequired = (function createBind2dTextureIfRequiredFunction(){
	var currentlyBoundTexture=null;
	return function(texToBind){	//todo specify texture index. maybe combine with setting texture index so can automatically keep textures loaded
								//curently just assuming using tex 0, already set as active texture (is set active texture a fast gl call?)
		if (texToBind != currentlyBoundTexture){
			gl.bindTexture(gl.TEXTURE_2D, texToBind);
			currentlyBoundTexture = texToBind;
		}
	}
})();

var mvMatrix = mat4.create();
//mat4.identity(mvMatrix);
var pMatrix = mat4.create();
mat4.identity(pMatrix);
//mat4.translate(mvMatrix,vec3.fromArray([0,0,-10])); //glmatix expects own vector type


var playerPosition = [-0.75,-0.5,-0.5];	//should be contained within playerMatrix. TODO sort this (code here consistent with webglwideanglecamera project)

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
	sparseSeparateFaces:false,
	sparseSmoothed:false,
	sparseWithNormals:false,	//TODO with/without color
	sparseBasicAvgSmoothed:false,
	sparseDCSmoothed:false,
	sparseDCWithNormals:false,	//TODO with/without color
	textured:false,
	normalMapped:true,
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

var logLimiter =0;
var globalAxisCollisionData;	//TODO less crappy way to return this.
var skipPolyToVoxStuff=true;

// do in n^2 checks if cast a series of parallel lines, assuming each starting from outside object. sort collisions, step through cells
var startGenFuncTime = Date.now();
var fromPolyModelFunctionFast = (function generateFromPolyModelFunctionFast(){
	
	if (skipPolyToVoxStuff){return;}	//TODO less crap way to skip this
	
	var selectedData;
	
	var teapotColData={};	
	var teapotObject = loadBlenderExport(teapotData);	//isn't actually a blender export - just a obj json
	loadColData(teapotColData, teapotObject);
	selectedData = teapotColData;
	
	/*
	var sshipColData={};	
	var sshipObject = loadBlenderExport(sshipData);	//isn't actually a blender export - just a obj json
	//shrink vertices
	var newVs = [];
	var oldVs = sshipObject.vertices;
	for (var vv in oldVs){
		newVs.push(oldVs[vv]/2.1);
	}
	sshipObject.vertices = newVs;
	loadColData(sshipColData, sshipObject);
	selectedData = sshipColData;
	*/
	
	var faces = selectedData.faces;
	var verts = selectedData.verts;
	
	//rotate teapot to aid debugging
	console.log("will attempt to rotate vert data");
	console.log(verts.length);
	console.log(verts[0]);
	var tmp;
	var thisVert;
	for (var vv=0;vv<verts.length;vv++){
		thisVert = verts[vv];
		tmp=thisVert[0];
		thisVert[0] = (thisVert[0]+thisVert[2])*0.7;
		thisVert[2] = (thisVert[2]-tmp)*0.7;
		thisVert[1] -=0.1;	//this line (moving teapot vertically) magically fixes things (wacky dual contouring* on half of teapot. i know not why. (*just basic intersection point averaging currently)
	}
	
	
	//precalculate face normal data
	/*
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
	*/
	//pregenerate and store data.
	//to use some assumptions used to generate smooth data, return bilinear smoothed data
	var myVoxData = [];
	var axisCollisionData={x:[],y:[],z:[]};	//data for collision point between grid points for some grid point to neighbouring point along given axis (where changes from voxel off/on)
	
	var myscale = 32;	//32 sort of works. other values don't. suppose because collisionPoints are *32 before returning
	
	for (var ii=0;ii<blocksize;ii++){
		var slicedata = [];
		myVoxData.push(slicedata);
		for (var jj=0;jj<blocksize;jj++){
			var stripdata = [];
			slicedata.push(stripdata);
			//var collisionData = checkForCollisions([(ii-30)/31,(jj-16)/31,-1],[(ii-30)/31,(jj-16)/31,1]);	//TODO rotate teapot 45 deg to fit into box better
			var collisionData = checkForCollisions([(ii-32)/myscale,(jj-32)/myscale,-1],[(ii-32)/myscale,(jj-32)/myscale,1]);	//TODO rotate teapot 45 deg to fit into box better
			collisionData.sort(function(a,b){return b.z-a.z;});	//sort collision data.
			var zidx=0;
			var nextCollision;
			var fill=-1;
			while (nextCollision = collisionData.pop()){
				if (nextCollision.z<0){console.log("nextCollision.z<0 - " + nextCollision.z);}
				axisCollisionData.z[ 64*64*ii + 64*jj + Math.floor(nextCollision.z) ] = nextCollision;	//could just keep nextCollision.z because x,y implied
				while (zidx<nextCollision.z){
					stripdata[zidx++]=fill;
				}
				fill = nextCollision.fill;	//transition to filled ie face downward or unfilled ie face up
			}
			while (zidx<blocksize){
				stripdata[zidx++]=fill;
			}
			
			//get collision data for x,y axes.
			//todo include normal data to do "proper" dual contouring instead of simple average of points
			collisionData = checkForCollisions([-1,(ii-32)/myscale,(jj-32)/myscale],[1,(ii-32)/myscale,(jj-32)/myscale]);
			//collisionData.sort(function(a,b){return b.x-a.x;});	//sort collision data.
			while (nextCollision = collisionData.pop()){
				if (nextCollision.x<0){console.log("nextCollision.x<0 - " + nextCollision.x);}
				axisCollisionData.x[ 64*64*Math.floor(nextCollision.x) + 64*ii + jj ] = nextCollision;	//could just keep nextCollision.x because y,z implied
			}
			collisionData = checkForCollisions([(ii-32)/myscale,-1,(jj-32)/myscale],[(ii-32)/myscale,1,(jj-32)/myscale]);
			//collisionData.sort(function(a,b){return b.y-a.y;});	//sort collision data.
			while (nextCollision = collisionData.pop()){
				if (nextCollision.y<0){console.log("nextCollision.y<0 - " + nextCollision.y);}
				axisCollisionData.y[ 64*64*ii + 64*Math.floor(nextCollision.y) + jj ] = nextCollision;	//could just keep nextCollision.y because x,z implied
			}
		}
	}
	console.log("myVoxData:");
	console.log(myVoxData);
	globalAxisCollisionData = axisCollisionData;
	console.log("axisCollisionData:");
	console.log(axisCollisionData);
	
	function checkForCollisions(raystart, rayend){

		//move teapot to see whether problems move with it
		//possibly should also shift output
		/*
		var indexToShift=1;
		var amountToShift=1/64;
		raystart[indexToShift]+=amountToShift;
		rayend[indexToShift]+=amountToShift;
	*/
		var collisionPoints = [];
		var rayVec=[ rayend[0]-raystart[0], rayend[1]-raystart[1], rayend[2]-raystart[2] ];
	
		for (var ff=0;ff<faces.length;ff++){
			var thisTri = faces[ff];
			
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
				
				/*
				//simplify above given that rayVec x,y components =0
				var displacement = [raystart[0]-thisVert[0],raystart[1]-thisVert[1]];
				var edgevector = [nextVert[0]-thisVert[0],nextVert[1]-thisVert[1]];
				var crossProdZ = displacement[0]*edgevector[1] - displacement[1]*edgevector[0];
				var dotProd = crossProdZ*rayVec[2];
				*/
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
				
				var triPoints = [verts[thisTri[0]], verts[thisTri[1]], verts[thisTri[2]]];
				var edgeVecs = [[triPoints[1][0] - triPoints[0][0], triPoints[1][1] - triPoints[0][1], triPoints[1][2] - triPoints[0][2]],
								[triPoints[2][0] - triPoints[0][0], triPoints[2][1] - triPoints[0][1], triPoints[2][2] - triPoints[0][2]]];
				var normVector = [edgeVecs[1][1]*edgeVecs[0][2] - edgeVecs[1][2]*edgeVecs[0][1],
									edgeVecs[1][2]*edgeVecs[0][0] - edgeVecs[1][0]*edgeVecs[0][2],
									edgeVecs[1][0]*edgeVecs[0][1] - edgeVecs[1][1]*edgeVecs[0][0]];
									
				
				//copied code from collision test project
				var vert0 = triPoints[0];
				var displacement = [raystart[0]-vert0[0],raystart[1]-vert0[1],raystart[2]-vert0[2]];
				var faceNormal = normVector;
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
								
									
					
				//collisionPoints.push({x:(colpoint[0]+1)*32, y:(colpoint[1]+1)*32, z:(colpoint[2]+1)*32, fill:signsSum<0?-1:1});	//TODO check polarity
				collisionPoints.push({x:(colpoint[0]+1)*32, y:(colpoint[1]+1)*32, z:(colpoint[2]+1)*32, fill:signsSum<0?-1:1});	//TODO check polarity
				
			}
		}
		
		if (collisionPoints.length>0 && logLimiter<10){
			console.log("created collision data. length = " + collisionPoints.length );
			console.log(collisionPoints);	//won't show in log because gets popped empty after returned from this func
			console.log(collisionPoints[0]);
			logLimiter++;
		}
		
		//return some dummy data
		//return [{fill:1,z:5},{fill:0,z:10}];
		return collisionPoints;
	}
	
	return function(ii,jj,kk){
		//avoid failure to run problem. TODO this properly
		ii = Math.max(0,ii);
		jj = Math.max(0,jj);
		kk = Math.max(0,kk);
		
		/*
		ii+=0.5;	//??
		jj+=0.5;	//??
		kk+=0.5;	//??
		*/
		
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
			if (typeof myVoxData == 'undefined'){console.log("myVoxData undefined");}
			if (typeof myVoxData[iif] == 'undefined'){console.log("myVoxData[" + iif + "] undefined");}
			if (typeof myVoxData[iif][jjf] == 'undefined'){console.log("myVoxData[" + iif + "][" + jjf +"] undefined");}
			if (typeof myVoxData[iif][jjf][kkf] == 'undefined'){console.log("myVoxData[" + iif + "][" + jjf +"][" + kkf + "] undefined");}
			*/
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
var voxFunction;

function init(){
	
	stats = new Stats();
	stats.showPanel( 0 ); // 0: fps, 1: ms, 2: mb, 3+: custom
	document.body.appendChild( stats.dom );
	
	var gui = new dat.GUI();
	gui.add(guiParams, "box");
	gui.add(guiParams, "stupid");
	gui.add(guiParams, "sparse");
	gui.add(guiParams, "sparseSeparateFaces");
	gui.add(guiParams, "sparseSmoothed");
	gui.add(guiParams, "sparseBasicAvgSmoothed");
	gui.add(guiParams, "sparseWithNormals");
	gui.add(guiParams, "sparseDCSmoothed");
	gui.add(guiParams, "sparseDCWithNormals");
	gui.add(guiParams, "textured");
	gui.add(guiParams, "normalMapped");
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
	
	//voxFunction = sinesfunction;
	voxFunction = sinesfunctiontwo;
	//voxFunction = landscapeFunction;
	//voxFunction = bigBallFunction;
	//voxFunction = bigCylinderFunction;
	//voxFunction = curveCornerFunction;
	//voxFunction = perlinPlanetFunction;
	//voxFunction = twistedTowerFunction;
	//voxFunction = roundedTwistedTowerFunction;
	//voxFunction = donutFunction;
	//voxFunction = bilinearFilterBinaryFunctionGen(sinesfunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(landscapeFunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(bigBallFunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(bigCylinderFunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(curveCornerFunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(perlinPlanetFunction);
	//voxFunction = bilinearFilterBinaryFunctionGen(twistedTowerFunction);
	
	seedValue= Math.random();
	//seedValue= 0.06941288120167277;	//interesting landscape with tunnels
	console.log("seed: " + seedValue);
	var genStartTime = Date.now();
	noise.seed(seedValue);
	//voxFunction = perlinfunction;
	//voxFunction = perlinfunctionTwoSided;
	//voxFunction = bilinearFilterBinaryFunctionGen(perlinfunction);
	
	
	//var voxFunction = fromPolyModelFunction;
	//var voxFunction = fromPolyModelFunctionFast;
	
	makeVoxdataForFunc(voxFunction);
	console.log("Time taken to generate: " + (Date.now()-genStartTime));
	
	function perlinfunction(ii,jj,kk){
		//return 10*noise.perlin3(ii/12,jj/12,kk/12);	//if divide by too small number, too many indices generated
		return 10*noise.perlin3(ii/24,jj/12,kk/12) - 0.75*kk +20;	//landscape with 3d perlin surface
	}
	function perlinfunctionTwoSided(ii,jj,kk){
		//return 10*noise.perlin3(ii/10,jj/10,kk/10) - 0.02*(kk-32)*(kk-32);	//landscape with 3d perlin surface
		return 10*wrapPerlin(ii/8,jj/8,kk/8,64/8) - 0.002*(kk-32)*(kk-32);	//landscape with 3d perlin surface
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
	function sinesfunctiontwo(ii,jj,kk){
		var sinscale=4/Math.PI;
		//return Math.sin(ii/sinscale)+Math.sin(jj/sinscale)+Math.sin(kk/sinscale);
		return Math.sin(ii/sinscale)+Math.sin(jj/sinscale)- kk/sinscale + 10;
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
	function twistedTowerFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		//kkm = kk-blocksize/2;
		var ang = kk*Math.PI/64;
		var iir = iim*Math.cos(ang)-jjm*Math.sin(ang);
		var jjr = iim*Math.sin(ang)+jjm*Math.cos(ang);
		//return 10 - Math.max(Math.abs(iir), Math.abs(jjr));
		var circleDist = Math.sqrt( Math.pow(Math.max(Math.abs(iir)-10 , 0),2) + Math.pow( Math.max(Math.abs(jjr)-10 , 0),2) );
		return circleDist>0? -circleDist : 10 - Math.max(Math.abs(iir), Math.abs(jjr));
	}
	function roundedTwistedTowerFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		//kkm = kk-blocksize/2;
		var ang = kk*Math.PI/64;
		var iir = iim*Math.cos(ang)-jjm*Math.sin(ang);
		var jjr = iim*Math.sin(ang)+jjm*Math.cos(ang);
		//return 10 - Math.max(Math.abs(iir), Math.abs(jjr));
		var circleDist = Math.sqrt( Math.pow(Math.max(Math.abs(iir)-10 , 0),2) + Math.pow( Math.max(Math.abs(jjr)-10 , 0),2) );
		//return circleDist>0? -circleDist : 10 - Math.max(Math.abs(iir), Math.abs(jjr));
		return circleDist>0? -circleDist+2 : 12 - Math.max(Math.abs(iir), Math.abs(jjr));
	}
	function donutFunction(ii,jj,kk){
		var iim,jjm,kkm;
		iim = ii-blocksize/2;
		jjm = jj-blocksize/2;
		kkm = kk-blocksize/2;
		var rad = Math.sqrt(iim*iim + jjm*jjm+0.001);
		var donutsradsq = (rad-21)*(rad-21) + kkm*kkm;
		return 100-donutsradsq;
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
	
	var mattoinvert = mat3.create();
	var myvec3 = vec3.create();
	
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
		var smoothVertices = [];
		var basicAvgVertices = [];
		var dcVertices = [];
		var normals = [];
		var dcNormals = [];
		var colors = [];
		var dcColors = [];
		var delta = 0.01;
		var badVertCount=0;
		
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
					if (ii>0){numNeighbours*=2;};
					if (jj>0){numNeighbours*=2;};
					if (kk>0){numNeighbours*=2;};
					
					if (voxdata[ii%64][jj%64][kk%64]){
						occNeighbours++;
					}
					if (kk>0 && voxdata[ii%64][jj%64][kk-1]){
						occNeighbours++;
					}
					if (jj>0 && voxdata[ii%64][jj-1][kk%64]){
						occNeighbours++;
					}
					if (jj>0 && kk>0 && voxdata[ii%64][jj-1][kk-1]){
						occNeighbours++;
					}
					if (ii>0 && voxdata[ii-1][jj%64][kk%64]){
						occNeighbours++;
					}
					if (ii>0 && kk>0 && voxdata[ii-1][jj%64][kk-1]){
						occNeighbours++;
					}
					if (ii>0 && jj>0 && voxdata[ii-1][jj-1][kk%64]){
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
			
			//info for dual contouring. TODO more efficient to put this with vertex creation, ie looking up voxdata[ii][jj][kk] etc...
			var ii_lo = ii-1;
			var jj_lo = jj-1;
			var kk_lo = kk-1;
			//TODO store evaluated function vals - already done this to determine whether each point is in/out
			var vfunc = voxFunction;	//to make more readable
			var vdata = [ vfunc(ii_lo,jj_lo,kk_lo), vfunc(ii_lo,jj_lo,kk), vfunc(ii_lo,jj,kk_lo), vfunc(ii_lo,jj,kk),
						vfunc(ii,jj_lo,kk_lo), vfunc(ii,jj_lo,kk), vfunc(ii,jj,kk_lo), vfunc(ii,jj,kk)];
						
			ii-=0.5;	//???
			jj-=0.5;
			kk-=0.5;
			
			
			vertices.push(ii/32, jj/32, kk/32);
			
			
			//normals.push(0,0,0);
				//^^ little faster for doing O(3) slow teapot thing

			//smooth vertices (TODO make 2 vert buffers and UI to switch between) 
			//just look at gradient between surrounding points, move downhill. or take analytic gradient from something function used to generate vox data
			//to make this work generally without needing to calculate analytic derivatives, just use numerical differences.
			
			var fudgeFactor = 0.5;	//less than 1 to avoid overshoot
			var centralPoint,gradX,gradY,gradZ,totalGradSq,sharedPart;
			
			for (var iter=0;iter<10;iter++){
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
			
			smoothVertices.push(ii/32, jj/32, kk/32);
			
			
			var invLengthSq = 1/( totalGradSq + 0.001)
			var invLength = Math.sqrt(invLengthSq);
			var normal = [invLength*gradX, invLength*gradY, invLength*gradZ];
			var normalOverLength = [invLengthSq*gradX, invLengthSq*gradY, invLengthSq*gradZ];
			normals.push(normal[0], normal[1], normal[2]);
			
			var grayColor = grayColorForPointAndNormal(ii,jj,kk,normal, invLength);	
			colors.push(grayColor, grayColor, grayColor);	//TODO separate surf color for directional lighting from ambient response
			
			
			//"dual contouring" ? 
			//AFAIK this is where find vertex positions by... find intersection of isosurface between neighbouring voxel centres of different polarity, for the 8 voxels around vertex position - 12 possible pairs. at these points, find surface normal. suppose that vertex position lies near to plane defined by this position and vertex. for gauss distribution about this point/plane, can multiply gaussian probabilities and find maximum/centre. equivalent to adding exponents. add some extra term to prefer points nerer centre of 8 nearby voxel centres (ie where vertex for simple minecraft style boxel would be).
			//p = point (should find p that minimises this func)
			//c = centre of given intersection point (plane centre)
			//n = plane normal
			//ie something like, minimise sum{( (p - c).n)^2 + const1*(p-c)*(p-c)} + const2*p*p (where origin at boxel vert pos)
			//don't need both const1 and const2 to be nonzero
			
			//for working see paper (todo scan/write up)
			// maximum where derivative is 0. for each component...
			
			//a step towards this - find each of 12 possible intersection points, average point and normals.
			//TODO note - maybe by storing unperturbed points, can use sharp shading across sharp creases?
			
			//TODO intersection points can be calculated up to 4 times. should calculate once, reuse.
			
			var sumx=0;
			var sumy=0;
			var sumz=0;
			
			var sumnx=0;
			var sumny=0;
			var sumnz=0;
			
			var sumnxx=0;
			var sumnxy=0;
			var sumnxz=0;
			var sumnyy=0;
			var sumnyz=0;
			var sumnzz=0;
			
			var centrebias = 0.001;	//k2 from paper working. play with this value. guess should scale with number of points averaged.
			var nk1 = 0.05;
			var sumnxx=centrebias + nk1;
			var sumnyy=centrebias + nk1;
			var sumnzz=centrebias + nk1;
			var sumnxy=0;
			var sumnxz=0;
			var sumnyz=0;
			
			var sumnxxcx=0;
			var sumnxycy=0;
			var sumnxzcz=0;
			var sumnyycy=0;
			var sumnxycx=0;
			var sumnyzcz=0;
			var sumnzzcz=0;
			var sumnxzcx=0;
			var sumnyzcy=0;
			
			var sumnum=0;
			//do sumx from lo,lo,lo point
			//todo iterative root funding
			//to calc normals (initially just check position finding
			//switch along z
			var intersect;
			var thisnorm;
			
			function addToSums(ii_rel,jj_rel,kk_rel){	//inputs are relative to lo corner
				if (isNaN(ii_rel) || isNaN(jj_rel) || isNaN(kk_rel)){
					console.log("nan input to addToSums!" + ii_rel + "," + jj_rel + "," + kk_rel);
				}
			
				var ii_fromcentre = ii_rel-0.5;
				var jj_fromcentre = jj_rel-0.5;
				var kk_fromcentre = kk_rel-0.5;
				var ii_world = ii_lo + ii_rel;
				var jj_world = jj_lo + jj_rel;
				var kk_world = kk_lo + kk_rel;
				
				//bodge?
				//ii_world+=0.5;
				//jj_world+=0.5;
				//kk_world+=0.5;
				
				//calculate normal for this position
				centralPoint = voxFunction(ii_world,jj_world,kk_world);
							
				gradX = (voxFunction(ii_world+delta,jj_world,kk_world)- centralPoint)/delta;
				gradY = (voxFunction(ii_world,jj_world+delta,kk_world)- centralPoint)/delta;
				gradZ = (voxFunction(ii_world,jj_world,kk_world+delta)- centralPoint)/delta;
				
				totalGradSq = gradX*gradX + gradY*gradY + gradZ*gradZ;
				
				//have some sort of hill normal. should move downhill
				//move by something like (discrepancy / totalGradent)*(gradientVector/totalGradient)
				//to avoid /0 error add something to totalGradient
			
				var totalgrad = Math.sqrt( totalGradSq + 0.000001);
				
				thisnorm = [ gradX/totalgrad , gradY/totalgrad , gradZ/totalgrad ];
				
				sumnx+=thisnorm[0];
				sumny+=thisnorm[1];
				sumnz+=thisnorm[2];
				
				sumnxx+=thisnorm[0]*thisnorm[0];
				sumnyy+=thisnorm[1]*thisnorm[1];
				sumnzz+=thisnorm[2]*thisnorm[2];
				sumnxy+=thisnorm[0]*thisnorm[1];
				sumnxz+=thisnorm[0]*thisnorm[2];
				sumnyz+=thisnorm[1]*thisnorm[2];
				
				sumnxxcx+=thisnorm[0]*thisnorm[0]*ii_fromcentre;	//TODO formulate using loop/matrix etc
				sumnxycy+=thisnorm[0]*thisnorm[1]*jj_fromcentre;
				sumnxzcz+=thisnorm[0]*thisnorm[2]*kk_fromcentre;
				sumnxycx+=thisnorm[0]*thisnorm[1]*ii_fromcentre;
				sumnyycy+=thisnorm[1]*thisnorm[1]*jj_fromcentre;
				sumnyzcz+=thisnorm[1]*thisnorm[2]*kk_fromcentre;
				sumnxzcx+=thisnorm[0]*thisnorm[2]*ii_fromcentre;
				sumnyzcy+=thisnorm[1]*thisnorm[2]*jj_fromcentre;
				sumnzzcz+=thisnorm[2]*thisnorm[2]*kk_fromcentre;
			}
			
			var valAtIntersect;
			
			if ( (vdata[0]-vdata[1]) && vdata[0]*vdata[1]<=0){	//sign switch from lo,lo,lo to lo,lo,hi
				sumnum++;
				intersect = vdata[0]/(vdata[0]-vdata[1]);
				intersect=modifyIntersect(intersect,ii_lo,jj_lo,kk_lo+intersect,0,1);
				sumz+= intersect;
				addToSums(0,0,intersect);
			}
			if ( (vdata[2]-vdata[3]) && vdata[2]*vdata[3]<=0){	//sign switch from lo,hi,lo to lo,hi,hi
				sumnum++;
				sumy++;
				intersect = vdata[2]/(vdata[2]-vdata[3]);
				intersect=modifyIntersect(intersect,ii_lo,jj_lo+1,kk_lo+intersect,2,3);
				sumz+= intersect;
				addToSums(0,1,intersect);
			}
			if ( (vdata[4]-vdata[5]) && vdata[4]*vdata[5]<=0){	//sign switch from hi,lo,lo to hi,lo,hi
				sumnum++;
				sumx++;
				intersect = vdata[4]/(vdata[4]-vdata[5]);
				intersect=modifyIntersect(intersect,ii_lo+1,jj_lo,kk_lo+intersect,4,5);
				sumz+= intersect;
				addToSums(1,0,intersect);
			}
			if ( (vdata[6]-vdata[7]) && vdata[6]*vdata[7]<=0){	//sign switch from hi,hi,lo to hi,hi,hi
				sumnum++;
				sumx++;
				sumy++;
				intersect = vdata[6]/(vdata[6]-vdata[7]);
				intersect=modifyIntersect(intersect,ii_lo+1,jj_lo+1,kk_lo+intersect,6,7);
				sumz+= intersect;
				addToSums(1,1,intersect);
			}
			//switch along y
			if ( (vdata[0]-vdata[2]) && vdata[0]*vdata[2]<=0){	//sign switch from lo,lo,lo to lo,hi,lo
				sumnum++;
				intersect = vdata[0]/(vdata[0]-vdata[2]);
				intersect=modifyIntersect(intersect,ii_lo,jj_lo+intersect,kk_lo,0,2);
				sumy+= intersect;
				addToSums(0,intersect,0);
			}
			if ( (vdata[1]-vdata[3]) && vdata[1]*vdata[3]<=0){	//sign switch from lo,lo,hi to lo,hi,hi
				sumnum++;
				sumz++;
				intersect = vdata[1]/(vdata[1]-vdata[3]);
				intersect=modifyIntersect(intersect,ii_lo,jj_lo+intersect,kk_lo+1,1,3);
				sumy+= intersect;
				addToSums(0,intersect,1);
			}
			if ( (vdata[4]-vdata[6]) && vdata[4]*vdata[6]<=0){	//sign switch from hi,lo,lo to hi,hi,lo
				sumnum++;
				sumx++;
				intersect =  vdata[4]/(vdata[4]-vdata[6]);
				intersect=modifyIntersect(intersect,ii_lo+1,jj_lo+intersect,kk_lo,4,6);
				sumy+= intersect;
				addToSums(1,intersect,0);
			}
			if ( (vdata[5]-vdata[7]) && vdata[5]*vdata[7]<=0){	//sign switch from hi,lo,hi to hi,hi,hi
				sumnum++;
				sumx++;
				sumz++;
				intersect =  vdata[5]/(vdata[5]-vdata[7]);
				intersect=modifyIntersect(intersect,ii_lo+1,jj_lo+intersect,kk_lo+1,5,7);
				sumy+= intersect;
				addToSums(1,intersect,1);
			}
			//switch along x
			if ( (vdata[0]-vdata[4]) && vdata[0]*vdata[4]<0){	//sign switch from lo,lo,lo to hi,lo,lo
				sumnum++;
				intersect = vdata[0]/(vdata[0]-vdata[4]);
				intersect=modifyIntersect(intersect,ii_lo+intersect,jj_lo,kk_lo,0,4);
				sumx+= intersect;
				addToSums(intersect,0,0);
			}
			if ( (vdata[1]-vdata[5]) && vdata[1]*vdata[5]<0){	//sign switch from lo,lo,hi to hi,lo,hi
				sumnum++;
				sumz++;
				intersect = vdata[1]/(vdata[1]-vdata[5]);
				intersect=modifyIntersect(intersect,ii_lo+intersect,jj_lo,kk_lo+1,1,5);
				sumx+= intersect;
				addToSums(intersect,0,1);
			}
			if ( (vdata[2]-vdata[6]) && vdata[2]*vdata[6]<0){	//sign switch from lo,hi,lo to hi,hi,lo
				sumnum++;
				sumy++;
				intersect = vdata[2]/(vdata[2]-vdata[6]);
				intersect=modifyIntersect(intersect,ii_lo+intersect,jj_lo+1,kk_lo,2,6);
				sumx+= intersect;
				addToSums(intersect,1,0);
			}
			if ( (vdata[3]-vdata[7]) && vdata[3]*vdata[7]<0){	//sign switch from lo,lo,hi to hi,lo,hi
				sumnum++;
				sumy++;
				sumz++;
				intersect = vdata[3]/(vdata[3]-vdata[7]);
				intersect=modifyIntersect(intersect,ii_lo+intersect,jj_lo+1,kk_lo+1,3,7);
				sumx+= intersect;
				addToSums(intersect,1,1);
			}
			
			function modifyIntersect(intersect,xx,yy,zz,idx_a,idx_b){	//try to smooth off artifacts along sharp edges. if not having any sharp edges, this is unnecessary perf drain
				//return intersect;	//turn off
				valAtIntersect = vfunc(xx,yy,zz);
				if (valAtIntersect==0){return intersect;}
				
				if (valAtIntersect*vdata[idx_a]<0){
					intersect*= vdata[idx_a]/(vdata[idx_a]-valAtIntersect);
					if (isNaN(intersect)){console.log("intersect is nan - if!!");}
				}else{
					intersect=intersect+(1-intersect)*valAtIntersect/(valAtIntersect-vdata[idx_b]);
				}
				if (isNaN(intersect)){console.log("intersect is nan!!");}
				return intersect;
			}
			
			//do matrix calculation using sums
			//glmatrix library provides a function to invert a matrix. can't find a func to multiply a vector by a matrix though, but simple to write func
				centrebias=0;	//to stop adding here
			mattoinvert[0]=sumnxx + centrebias*sumnum;	mattoinvert[1]=sumnxy;						mattoinvert[2]=sumnxz;	
			mattoinvert[3]=sumnxy;						mattoinvert[4]=sumnyy + centrebias*sumnum;	mattoinvert[5]=sumnyz;
			mattoinvert[6]=sumnxz;						mattoinvert[7]=sumnyz;						mattoinvert[8]=sumnzz + centrebias*sumnum;
			
			for (var cc=0;cc<9;cc++){
				if (isNaN(mattoinvert[cc])){
					console.log("NaN!!" + mattoinvert);
				}
			}
			
			mattoinvert = mat3.inverse(mattoinvert);	//since matrix is symmetric, can inversion be done more efficiently?
			
			myvec3[0] = sumnxxcx + sumnxycy + sumnxzcz + nk1*(sumx/sumnum - 0.5);
			myvec3[1] = sumnxycx + sumnyycy + sumnyzcz + nk1*(sumy/sumnum - 0.5);
			myvec3[2] = sumnxzcx + sumnyzcy + sumnzzcz + nk1*(sumz/sumnum - 0.5);
			
			if (mattoinvert==null){console.log("null matrix!!!!");};	//else{console.log("matrix is not null");}
			mat3.multiplyVec3(mattoinvert, myvec3);
			
			var dcPos = [ii_lo+myvec3[0], jj_lo+myvec3[1], kk_lo+myvec3[2]];	
			dcVertices.push(dcPos[0]/32, dcPos[1]/32, dcPos[2]/32);
			
			//normalise average normals
			var sumNorm = Math.sqrt(sumnx*sumnx + sumny*sumny + sumnz*sumnz + 0.1);
			var dcNorm = [sumnx/sumNorm, sumny/sumNorm, sumnz/sumNorm];
			dcNormals.push( dcNorm[0], dcNorm[1], dcNorm[2]);
			
			/*
			//basic average points method
			if (sumnum>0){	//AFAIK this should not happen
				ii_lo+=sumx/sumnum;
				jj_lo+=sumy/sumnum;
				kk_lo+=sumz/sumnum;
			}
			*/
			var sums;
			
			//override basic average vertex if have intersection data from poly->vox (todo skip preceding code in this case)
			if (typeof globalAxisCollisionData =="undefined"){
				basicAvgVertices.push(ii_lo/32, jj_lo/32, kk_lo/32);
			}else{
				
				var ii_adj = ii_lo+0;	//adjust to guess fix
				var jj_adj = jj_lo+0;
				var kk_adj = kk_lo+0;
				
				var idxToLookup = 64*64*ii_adj + 64*jj_adj + kk_adj;			//todo check out by ones
				sums = {x:0,y:0,z:0,n:0};
	
				addCollisionData(globalAxisCollisionData.x[idxToLookup]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup + 64]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup + 1]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup + 64+1]);
				
				addCollisionData(globalAxisCollisionData.y[idxToLookup]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup + 64*64]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup + 1]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup + 64*64+1]);
				
				addCollisionData(globalAxisCollisionData.z[idxToLookup]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup + 64*64]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup + 64]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup + 64*64+64]);
				
				/*
				//guess at correct?
				addCollisionData(globalAxisCollisionData.x[idxToLookup]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup - 64]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup - 1]);
				addCollisionData(globalAxisCollisionData.x[idxToLookup - 64+1]);
				
				addCollisionData(globalAxisCollisionData.y[idxToLookup]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup - 64*64]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup - 1]);
				addCollisionData(globalAxisCollisionData.y[idxToLookup - 64*64+1]);
				
				
				addCollisionData(globalAxisCollisionData.z[idxToLookup]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup - 64*64]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup - 64]);
				addCollisionData(globalAxisCollisionData.z[idxToLookup - 64*64+64]);
				*/
				
				if (sums.n ==0){
					alert("sums.n = 0 !!");
					//basicAvgVertices.push(0, 0, 0);	//THIS SHOULD NOT HAPPEN!!!
					basicAvgVertices.push(1, 1, 1);	//THIS SHOULD NOT HAPPEN!!!
					badVertCount++;
				}
				//console.log(sums);
				//basicAvgVertices.push(sums.x/sums.n, sums.y/sums.n, sums.z/sums.n);
				basicAvgVertices.push((sums.x/sums.n)/32, (sums.y/sums.n)/32, (sums.z/sums.n)/32);
			}
				//TODO move function outside loop!!!
			function addCollisionData(cdata){
				if (typeof cdata == "undefined"){return;}
				sums.x+= cdata.x;
				sums.y+= cdata.y;
				sums.z+= cdata.z;
				sums.n++;
			}
			
			
			
			
			//grayColor = grayColorForPointAndNormal(dcPos[0],dcPos[1],dcPos[2],dcNorm,1/sumNorm);	//wierd result since average normal is not the normal at this point! 
			grayColor = grayColorForPointAndNormal(dcPos[0],dcPos[1],dcPos[2]);	//don't pass in normal info, calc inside function
			
			dcColors.push(grayColor, grayColor, grayColor);	//TODO separate surf color for directional lighting from ambient response
		}

		
		
		function grayColorForPointAndNormal(ii,jj,kk, normal, invLength){	//note can just calculate normal at point, but saves some calculation if already have it
		
			//guess
			var fudgeOffset = 0.5;	//TODO work out correct value!
			
			ii+=fudgeOffset;
			jj+=fudgeOffset;
			kk+=fudgeOffset;
		
			//colors.push(Math.random(),Math.random(),Math.random());
			var colorScale = 6;	//scale of noise (smaller = finer)
			//var grayColor = 0.5+0.5*noise.perlin3(ii/colorScale,jj/colorScale,kk/colorScale);	// mapt -1 to 1 -> 0 to 1
			//var grayColor = 1.5+0.5*noise.perlin3(ii/colorScale,jj/colorScale,kk/colorScale);
			var grayColor = 1+1*sumPerlin(ii/colorScale,jj/colorScale,kk/colorScale,1.8);
			grayColor*=0.5;
			//colors.push(grayColor, grayColor, grayColor);
			//normals.push(0,0,0); 
			
			//colour by local curvature. guess an equation for this.
			//really a saddle should be more shady than a flat plane, and direction of lighting could be used, but simple version may provide ok effect 
			var twiceCentralPoint = 2*voxFunction(ii,jj,kk);
			var fwdX = voxFunction(ii+delta,jj,kk);
			var bwdX = voxFunction(ii-delta,jj,kk);
			var fwdY = voxFunction(ii,jj+delta,kk);
			var bwdY = voxFunction(ii,jj-delta,kk);
			var fwdZ = voxFunction(ii,jj,kk+delta);
			var bwdZ = voxFunction(ii,jj,kk-delta);
			
			var shiftX = fwdX + bwdX - twiceCentralPoint;
			var shiftY = fwdY + bwdY - twiceCentralPoint;
			var shiftZ = fwdZ + bwdZ - twiceCentralPoint;
			
			if (!normal){	//override input normal/invLength (if works, remove inputs)
				normal = [(fwdX-bwdX)/delta, (fwdY-bwdY)/delta, (fwdZ-bwdZ)/delta];	//over delta avoids numerical error afaik
				invLength = 1/Math.sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2] + 0.00001);
				normal = normal.map(function(elem){return elem*invLength;});	
			}
			
			//try laplacian. suspect should use with second derivative in direction of normal
			var curveColor = shiftX + shiftY + shiftZ;
			var shiftGrad = voxFunction(ii+delta*normal[0],jj+delta*normal[1],kk+delta*normal[2]) + voxFunction(ii-delta*normal[0],jj-delta*normal[1],kk-delta*normal[2]) - twiceCentralPoint;	//TODO formulate this using the existing samples?
			curveColor -= shiftGrad;
			
			//positive curvature doesn't increase lighting (unless saddle-like). TODO nonlinear shading(curve)
			curveColor = Math.max(curveColor,0.0);
			
			//divide by steepness of scalar field func
			curveColor*=invLength;
			curveColor/=delta*delta;
			
			curveColor = Math.max(Math.atan(curveColor)*(2/Math.PI),0.0);
			
			//grayColor*=0.5-curveColor;	//using *= to retain perlin
			grayColor=0.5-curveColor;	//using *= to retain perlin
			return grayColor;
		}

		function getNumberOfGridPoint(ii,jj,kk){
			return ii*65*65 + jj*65 + kk;
		}
		function pushIndexForNumber(indexarr,nii,njj,nkk){
			indexarr.push(indexForGridPoint[nii],indexForGridPoint[njj],indexForGridPoint[nkk]);
		}
		
		var oneVertIdx;
		var directionalIndices = [[],[],[],[],[],[]];
		for (var ii=0;ii<64;ii++){
			for (var jj=0;jj<64;jj++){
				for (var kk=0;kk<64;kk++){
					var difference=voxdata[ii][jj][(kk+1)%64]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = getNumberOfGridPoint(ii,jj,kk+1);
						if ( difference>0 ){
							pushIndexForNumber(directionalIndices[0], oneVertIdx , oneVertIdx+65 , oneVertIdx+65+65*65 );	//bottom faces
							pushIndexForNumber(directionalIndices[0],  oneVertIdx , oneVertIdx+65+65*65, oneVertIdx+65*65);
						}else{
							pushIndexForNumber(directionalIndices[1], oneVertIdx, oneVertIdx+65+65*65, oneVertIdx+65 );		//top faces
							pushIndexForNumber(directionalIndices[1], oneVertIdx, oneVertIdx+65*65, oneVertIdx+65+65*65);
							topFaceCount++;
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<64;ii++){
			for (var jj=0;jj<64;jj++){
				for (var kk=0;kk<64;kk++){
					var difference=voxdata[ii][(jj+1)%64][kk]-voxdata[ii][jj][kk];
					if (difference!=0){
						oneVertIdx = getNumberOfGridPoint(ii,jj+1,kk);
						if ( difference<0 ){
							pushIndexForNumber(directionalIndices[2], oneVertIdx, oneVertIdx+1, oneVertIdx+1+65*65 );
							pushIndexForNumber(directionalIndices[2], oneVertIdx, oneVertIdx+1+65*65, oneVertIdx+65*65);
						}else{
							pushIndexForNumber(directionalIndices[3], oneVertIdx, oneVertIdx+1+65*65, oneVertIdx+1 );
							pushIndexForNumber(directionalIndices[3], oneVertIdx, oneVertIdx+65*65, oneVertIdx+1+65*65);
						}
					}
				}
			}
		}
		
		for (var ii=0;ii<64;ii++){
			for (var jj=0;jj<64;jj++){
				for (var kk=0;kk<64;kk++){
					var difference=voxdata[(ii+1)%64][jj][kk]-voxdata[ii][jj][kk];
					if ( difference!=0 ){
						oneVertIdx = getNumberOfGridPoint(ii+1,jj,kk);					
						if ( difference>0 ){
							pushIndexForNumber(directionalIndices[4], oneVertIdx, oneVertIdx+1, oneVertIdx+1+65);
							pushIndexForNumber(directionalIndices[4], oneVertIdx, oneVertIdx+1+65, oneVertIdx+65);
						}else{
							pushIndexForNumber(directionalIndices[5], oneVertIdx, oneVertIdx+1+65, oneVertIdx+1);
							pushIndexForNumber(directionalIndices[5], oneVertIdx, oneVertIdx+65, oneVertIdx+1+65);
						}
					}
				}
			}
		}
		
		console.log("indices info:");
		console.log(directionalIndices[0].length);
		console.log(directionalIndices[1].length);
		console.log(directionalIndices[2].length);
		console.log(directionalIndices[3].length);
		console.log(directionalIndices[4].length);
		console.log(directionalIndices[5].length);
		
		console.log("badVertCount: " + badVertCount);
		
		return {
			vertices:vertices,
			smoothVertices:smoothVertices,
			basicAvgVertices:basicAvgVertices,
			dcVertices:dcVertices,
			normals:normals,
			dcNormals:dcNormals,
			colors:colors,
			dcColors:dcColors,
			directionalIndices:directionalIndices,
			indices:Array.prototype.concat.apply([],directionalIndices)
		};
	})();
	
	
	console.log("sparse verts : " + voxData["sparse"].vertices.length );
	
	console.log("topFaceCount: " + topFaceCount);
	
	
	initGL();
	
	gl.clearColor.apply(gl,[0.7,1,1,1]);
	
	initShaders();
	initTexture();
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


function drawObjectFromPreppedBuffers(bufferObj, shaderProg, startIdx, numIndices){
	gl.uniformMatrix4fv(shaderProg.uniforms.uPMatrix, false, pMatrix);	//TODO not every frame.
	gl.uniformMatrix4fv(shaderProg.uniforms.uMVMatrix, false, mvMatrix);
	
	if (typeof startIdx == "undefined"){
		gl.drawElements(gl.TRIANGLES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
		//gl.drawElements(gl.LINES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
	}else{
		gl.drawElements(gl.TRIANGLES, numIndices, gl.UNSIGNED_SHORT, 2*startIdx);	//unsure where the 2 comes from!
	}
	
	//make a 2x2 grid
	mat4.translate(mvMatrix, [0,2,0]);
	gl.uniformMatrix4fv(shaderProg.uniforms.uMVMatrix, false, mvMatrix);
	gl.drawElements(gl.TRIANGLES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
	mat4.translate(mvMatrix, [2,-2,0]);
	gl.uniformMatrix4fv(shaderProg.uniforms.uMVMatrix, false, mvMatrix);
	gl.drawElements(gl.TRIANGLES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
	mat4.translate(mvMatrix, [0,2,0]);
	gl.uniformMatrix4fv(shaderProg.uniforms.uMVMatrix, false, mvMatrix);
	gl.drawElements(gl.TRIANGLES, bufferObj.vertexIndexBuffer.numItems, gl.UNSIGNED_SHORT, 0);
	mat4.translate(mvMatrix, [-2,-2,0]);
}

var currentTime=0;

//borrowed from rayMarcherTest
var camSpeed = [0,0,0];	//TODO have collision test inside iterate mechanics
var iterateMechanics = (function generateIterateMechanics(){
	var newTime = Date.now();
	
	var wasdPressed=[false,false,false,false];
	
	window.addEventListener("keyup",function(e){
		var keyCode=e.keyCode;
		if (keyCode == 87){wasdPressed[0]=false;}
		if (keyCode == 65){wasdPressed[1]=false;}
		if (keyCode == 83){wasdPressed[2]=false;}
		if (keyCode == 68){wasdPressed[3]=false;}
	});
	
	window.addEventListener("keydown",function(e){
		var keyCode=e.keyCode;
		//console.log(keyCode);
		//WASD = 87,65,83,68
		
		if (keyCode == 87){wasdPressed[0]=true;}
		if (keyCode == 65){wasdPressed[1]=true;}
		if (keyCode == 83){wasdPressed[2]=true;}
		if (keyCode == 68){wasdPressed[3]=true;}
		
		var rollSpd = 0.05;
		var rollAmt = 0;
		if (keyCode == 69){rollAmt-=rollSpd;}
		if (keyCode == 81){rollAmt+=rollSpd;}
		rotatePlayer([0,0,-rollAmt]);	// - because suspect matrix inverted compared to rayMarcherTest...
	});
	
	//todo set keypresses to false on lose focus (exit page)
	
	var timeRemainder = 0;
	var camStep = 0.0002;
	
	return function(){
		var oldTime = newTime;
		newTime = Date.now();
		var timeElapsed = Math.min(newTime - oldTime, 1000);	//1s max
		
		//cap time elapsed
		timeElapsed = Math.min(timeElapsed,100);
		timeRemainder +=timeElapsed;
		
		while (timeRemainder > 0){
	
			var camMove = [0,0,0];
			if (wasdPressed[0]){camMove[2]+=camStep;}
			if (wasdPressed[1]){camMove[0]+=camStep;}
			if (wasdPressed[2]){camMove[2]-=camStep;}
			if (wasdPressed[3]){camMove[0]-=camStep;}
			
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
			
			camSpeed = camSpeed.map(function(component){return component*0.95;});	//linear drag
			
			playerPosition = playerPosition.map(function(val,ii){return val+camSpeed[ii];});
			//playerPosition = playerPosition.map(function(val,ii){return val;});
			
			timeRemainder-=10;
		}
		
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
	
	//document.getElementById("debugtext").innerHTML = checkCollision(playerPosition);
	document.getElementById("debugtext").innerHTML = checkCollisionForVoxFunction(playerPosition, 1);
}

var cubeSideLighting=[
	[1,0,0],
	[0,1,0],
	[0,0,1],
	[0,1,1],
	[1,0,1],
	[1,1,0]
];	//test colours. todo dot product face direction with light, make consistent with normals lighting
	

function drawScene(drawTime){
	//requestAnimationFrame(drawScene);
	
	stats.end();
	stats.begin();
	
	resizecanvas(1);	//TODO should this really happen every frame? perf impact?
	gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);	//TODO should this be every frame?

	mat4.perspective(mainCamFov, gl.viewportWidth/gl.viewportHeight, 0.001,20.0,pMatrix);	//apparently 0.9.5, last param is matrix rather than 1st!! todo use newer???
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
	if (guiParams.sparseSmoothed){
		drawObjectFromBuffers(voxBuffers["sparseSmoothed"], activeShaderProgram);
	}
	if (guiParams.sparseBasicAvgSmoothed){
		drawObjectFromBuffers(voxBuffers["sparseBasicAvgSmoothed"], activeShaderProgram);
	}
	if (guiParams.sparseDCSmoothed){
		drawObjectFromBuffers(voxBuffers["sparseDCSmoothed"], activeShaderProgram);
	}
	
	if (guiParams.sparseSeparateFaces){
		activeShaderProgram = shaderPrograms.simpleColor;
		gl.useProgram(activeShaderProgram);
		prepBuffersForDrawing(voxBuffers["sparse"], activeShaderProgram);
		var runningIndex=0;
		var sparseBuffers = voxBuffers["sparse"];
		for (var dd=0;dd<6;dd++){
			gl.uniform3fv(activeShaderProgram.uniforms.uColor, cubeSideLighting[dd]);
			drawObjectFromPreppedBuffers(sparseBuffers, activeShaderProgram, runningIndex, sparseBuffers.directionalIndices[dd].length);	//TODO just store lengths
			runningIndex+=sparseBuffers.directionalIndices[dd].length;
		}
	}
	
	if (guiParams.sparseWithNormals){
		//activeShaderProgram = shaderPrograms.withNormals;
		activeShaderProgram = shaderPrograms.withNormalsAndColor;
		gl.useProgram(activeShaderProgram);
		drawObjectFromBuffers(voxBuffers["sparseWithNormals"], activeShaderProgram);
	}
	if (guiParams.sparseDCWithNormals){
		//activeShaderProgram = shaderPrograms.withNormals;
		activeShaderProgram = shaderPrograms.withNormalsAndColor;
		gl.useProgram(activeShaderProgram);
		drawObjectFromBuffers(voxBuffers["sparseDCWithNormals"], activeShaderProgram);
	}
	
	if (guiParams.textured){
		activeShaderProgram = shaderPrograms.texmap;
		gl.useProgram(activeShaderProgram);
		bind2dTextureIfRequired(texture);
		gl.uniform1i(activeShaderProgram.uniforms.uSampler, 0);
		drawObjectFromBuffers(voxBuffers["sparseDCWithNormals"], activeShaderProgram);
	}
	
	if (guiParams.normalMapped){
		activeShaderProgram = shaderPrograms.nmap;
		gl.useProgram(activeShaderProgram);
		bind2dTextureIfRequired(texture);
		gl.uniform1i(activeShaderProgram.uniforms.uSampler, 0);
		drawObjectFromBuffers(voxBuffers["sparseDCWithNormals"], activeShaderProgram);
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
function checkCollisionForVoxFunction(position, size){
	//return voxFunction(-32*position[0],-32*position[1],-32*position[2])>0;
	var xx = -32*position[0];
	var yy = -32*position[1];
	var zz = -32*position[2];
	var delta = 0.01;
	var centreVal = voxFunction(xx,yy,zz);
	var gradX = voxFunction(xx+delta,yy,zz) - centreVal;
	var gradY = voxFunction(xx,yy+delta,zz) - centreVal;
	var gradZ = voxFunction(xx,yy,zz+delta) - centreVal;
	
	var gradLengthsq = gradX*gradX + gradY*gradY + gradZ*gradZ;
	var collisionPenetration = size/delta+centreVal/Math.sqrt(gradLengthsq);
	if (collisionPenetration>0){
		//reflect velocity. dot gradient with speed to check travelling towards surface.
		var gradDotSpeed = camSpeed[0]*gradX + camSpeed[1]*gradY + camSpeed[2]*gradZ;
		if ( gradDotSpeed<0){
			//new vel = vel + 2* normal.vel * normal
			var multiplier = 2*gradDotSpeed/gradLengthsq;
			camSpeed[0]-= multiplier*gradX;
			camSpeed[1]-= multiplier*gradY;
			camSpeed[2]-= multiplier*gradZ;
		}
	}
	
	return collisionPenetration>0;
}

function sumPerlin(ii,jj,kk,amplscale){
	//seems perlin lib doesn't have many options. want something with discernable texture over more length scales.
	var total=0;
	var colorScale=6;
	var amplitude=1.5;
	for (var iter=0;iter<10;iter++){
		colorScale/=2;
		amplitude/=amplscale;	//TODO sum series to normalise sum of amplitudes
		total+=amplitude*noise.perlin3(ii/colorScale,jj/colorScale,kk/colorScale);	//TODO consistent random offsets for levels (so doesn't spike at 0)
	}
	return total;
}

//seems like perlin library using does not wrap. TODO for cleanliness write own perlin (wrapping should be fairly easy) 
//for now bodge averaging 8 samples
function wrapPerlin(ii,jj,kk,wrapscale){
	ii%=wrapscale;	//this doesn't handle negative values
	jj%=wrapscale;
	kk%=wrapscale;
	
	var ii_fract = ii/wrapscale;
	var jj_fract = jj/wrapscale;
	var kk_fract = kk/wrapscale;
	
	var ii_otherfract = 1-ii_fract;
	var jj_otherfract = 1-jj_fract;
	var kk_otherfract = 1-kk_fract;
	
	var sum=0;
	
	sum+=noise.perlin3(ii, jj, kk)* ii_fract*jj_fract*kk_fract;
	sum+=noise.perlin3(ii, jj, kk+wrapscale)* ii_fract*jj_fract*kk_otherfract;
	sum+=noise.perlin3(ii, jj+wrapscale, kk)* ii_fract*jj_otherfract*kk_fract;
	sum+=noise.perlin3(ii, jj+wrapscale, kk+wrapscale)* ii_fract*jj_otherfract*kk_otherfract;
	sum+=noise.perlin3(ii+wrapscale, jj, kk)* ii_otherfract*jj_fract*kk_fract;
	sum+=noise.perlin3(ii+wrapscale, jj, kk+wrapscale)* ii_otherfract*jj_fract*kk_otherfract;
	sum+=noise.perlin3(ii+wrapscale, jj+wrapscale, kk)* ii_otherfract*jj_otherfract*kk_fract;
	sum+=noise.perlin3(ii+wrapscale, jj+wrapscale, kk+wrapscale)* ii_otherfract*jj_otherfract*kk_otherfract;
	
	return sum/8;	//not equivalent to wrapping perlin - generally result will be smaller.
}