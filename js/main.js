
var canvas;
//var gl;
var crosssectionscanvas;
var ctx;
var voxdata;


var basicCubeData = {
	vertices : [
		0,0,0,	//0
		0,0,1,	//1
		0,1,0,	//2
		0,1,1,	//3
		1,0,0,	//4
		1,0,1,	//5
		1,1,0,	//6
		1,1,1],	//7
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
}

var mvMatrix = mat4.create();
mat4.identity(mvMatrix);
var pMatrix = mat4.create();
mat4.identity(pMatrix);
//mat4.translate(mvMatrix,vec3.fromArray([0,0,-10])); //glmatix expects own vector type

mvMatrix[14]=-3;	//move back to look at thing (is this camera or thing position?)
mvMatrix[13]=-0.5;
mvMatrix[12]=0;

mat4.rotateX(mvMatrix, -1.3);	//rads
mat4.rotateZ(mvMatrix, 4);

voxdata = [];

var guiParams={
	box:false,
	stupid:false,
	sparse:true,
	sparseWithNormals:false
};

function init(){
	stats = new Stats();
	stats.showPanel( 0 ); // 0: fps, 1: ms, 2: mb, 3+: custom
	document.body.appendChild( stats.dom );
	
	var gui = new dat.GUI();
	gui.add(guiParams, "box");
	gui.add(guiParams, "stupid");
	gui.add(guiParams, "sparse");
	gui.add(guiParams, "sparseWithNormals");
	
	canvas=document.getElementById("glcanvas");
	gl=glcanvas.getContext("webgl");
	
	crosssectionscanvas=document.getElementById("mycanvas");
	ctx=crosssectionscanvas.getContext("2d");
	
	
	//voxel data. maybe better to have octree or similar. for now, use array of arrays. implementation maybe sparse. maybe good to use bits in a number (32 bit?)
	var blocksize = 64;	//64x64x64 = 512x512, which seems manageable to draw slices onto canvas
	var canvassize = 512;
	
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
	
	var voxFunction = sinesfunction;
	//var voxFunction = landscapeFunction;
	//var voxFunction = bigBallFunction;
	
	noise.seed(Math.random());
	//var voxFunction = perlinfunction;
	
	makeVoxdataForFunc(voxFunction);	
	
	function perlinfunction(ii,jj,kk){
		//return 10*noise.perlin3(ii/12,jj/12,kk/12);	//if divide by too small number, too many indices generated
		return 10*noise.perlin3(ii/24,jj/12,kk/12) - 0.75*kk +20;	//landscape with 3d perlin surface
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
			//vertices.push(ii/32, jj/32, kk/32);
			
			//smooth vertices (TODO make 2 vert buffers and UI to switch between) 
			//just look at gradient between surrounding points, move downhill. or take analytic gradient from something function used to generate vox data
			//to make this work generally without needing to calculate analytic derivatives, just use numerical differences.
			var delta = 0.01;
			
			var centralPoint = voxFunction(ii,jj,kk);
						
			var gradX = (voxFunction(ii+delta,jj,kk)- centralPoint)/delta;
			var gradY = (voxFunction(ii,jj+delta,kk)- centralPoint)/delta;
			var gradZ = (voxFunction(ii,jj,kk+delta)- centralPoint)/delta;
			
			var totalGradSq = gradX*gradX + gradY*gradY + gradZ*gradZ;
			
			//have some sort of hill normal. should move downhill
			//move by something like (discrepancy / totalGradent)*(gradientVector/totalGradient)
			//to avoid /0 error add something to totalGradient
		
			var sharedPart = centralPoint / ( totalGradSq + 0.01);
			
			var newi = ii-sharedPart*gradX;
			var newj = jj-sharedPart*gradY;
			var newk = kk-sharedPart*gradZ;
			
			vertices.push(newi/32, newj/32, newk/32);

			
			//do another normal measurement at the displaced point
			centralPoint = voxFunction(newi,newj,newk);
			gradX = (voxFunction(newi+delta,newj,newk)- centralPoint)/delta;
			gradY = (voxFunction(newi,newj+delta,newk)- centralPoint)/delta;
			gradZ = (voxFunction(newi,newj,newk+delta)- centralPoint)/delta;
			totalGradSq = gradX*gradX + gradY*gradY + gradZ*gradZ;
			
			var invLength = Math.sqrt(1/( totalGradSq + 0.001));
			normals.push(invLength*gradX, invLength*gradY, invLength*gradZ); 
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
	gl.useProgram(shaderPrograms.withNormals);
	drawObjectFromBuffers(voxBuffers["sparseWithNormals"], shaderPrograms.withNormals);
	
	gl.enable(gl.DEPTH_TEST);
	
	requestAnimationFrame(drawScene);
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
function drawScene(drawTime){
	requestAnimationFrame(drawScene);
	
	stats.end();
	stats.begin();
	
	resizecanvas(1);	//TODO should this really happen every frame? perf impact?
	gl.viewport(0, 0, gl.viewportWidth, gl.viewportHeight);	//TODO should this be every frame?

	mat4.perspective(60, gl.viewportWidth/gl.viewportHeight, 0.1,200.0,pMatrix);	//apparently 0.9.5, last param is matrix rather than 1st!! todo use newer???
																	//also old one uses degs!
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
		activeShaderProgram = shaderPrograms.withNormals;
		gl.useProgram(activeShaderProgram);
		drawObjectFromBuffers(voxBuffers["sparseWithNormals"], activeShaderProgram);
	}
	
	mat4.rotateZ(mvMatrix,0.0003*(drawTime-currentTime)); 
	currentTime=drawTime;
}
