window.addEventListener("load", main);

var results;
var canvasRef;
var ctxRef;
var limit = 1150;

function draw() {
   // Get the canvas reference
   canvasRef = document.getElementById('canvasRef');
   ctxRef = canvasRef.getContext('2d');
   ctxRef.clearRect(0, 0, canvasRef.width, canvasRef.height);

   var yOffset = 10; // padding between lines
   var yPos = 5; // Initial vertical position

   for (var key in results) {
      var obj = results[key];

      if (obj.hasOwnProperty('occurences') && Array.isArray(obj.occurences)) {
         
         // Actual drawing
         for (var i = 0; i < obj.occurences.length; i++) {
            if (obj.occurences[i].hasOwnProperty('position')) {
               var x = limit + obj.occurences[i].position; // x position
               var y = yPos + (i * yOffset); // y position
               
               // Draws occurrences
               ctxRef.fillStyle = 'blue';
               ctxRef.fillRect(x, y, 10, 5);

               // Draws sequence name
               ctxRef.fillStyle = 'black';
               ctxRef.fillText(obj.occurences[i].arn, 80, y + 5);

               // Draws the sequence line
               ctxRef.strokeStyle = 'black';
               ctxRef.beginPath();
               ctxRef.moveTo(150, y + 5);
               ctxRef.lineTo(limit, y + 5);
               ctxRef.lineTo(limit, y);
               ctxRef.lineTo(limit + 5 , y);
               ctxRef.stroke();
            }
         }
      }
      
      yPos += (obj.occurences.length + 1) * yOffset; // Adjusts the vertical position for the next line

      // Calculate the window size for drawing it
      windowStart = limit + obj.debut;
      windowEnd = limit + obj.fin + 1;
      windowSize = windowEnd - windowStart;

      // Draws the window
      var windowY = yPos-5;
      var windowHeight = yPos - 10;
      ctxRef.strokeStyle = 'red';
      ctxRef.lineWidth = 1;
      ctxRef.strokeRect(windowStart, windowY, windowEnd - windowStart, windowHeight);
   }

   // Function to display the occurence data
   canvasRef.addEventListener('click', function(event) {
      var rect = canvasRef.getBoundingClientRect();
      var mouseX = event.clientX - rect.left;
      var mouseY = event.clientY - rect.top;
   
      var yPos = 5;
   
      //Verify if the click is on the rectangle
      for (var key in results) {
         var obj = results[key];
         if (obj.hasOwnProperty('occurences') && Array.isArray(obj.occurences)) {
            for (var i = 0; i < obj.occurences.length; i++) {
               var x = limit + obj.occurences[i].position;
               var y = yPos + (i * yOffset);
   
               var rectWidth = 10;
               var rectHeight = 5;
               if (mouseX >= x && mouseX <= x + rectWidth && mouseY >= y && mouseY <= y + rectHeight) {
                  // Display info
                  alert('Sequence : ' + obj.occurences[i].arn +
                        '\nPosition : ' + obj.occurences[i].position +
                        '\nScore : '+ obj.occurences[i].score);
                  return;
               }
            }
         }
         yPos += (obj.occurences.length + 1) * yOffset;
      }
   });
   
}

function main() {
   //Selection du fichier Json 
   document.getElementById('inputfile')
      .addEventListener('change', function () {
         var fr = new FileReader();
         fr.onload = function () {
            results = JSON.parse(fr.result);
            draw();
         }
         fr.readAsText(this.files[0]);
      })
}
