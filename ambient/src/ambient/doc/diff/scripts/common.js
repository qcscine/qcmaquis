keys_height = 0;

function adjustScreen(){
    var height = Math.min(($(window).height()), parseInt($("section#intro").css("max-height"),10));
    $("section#intro").css({height: height});
    $("div#cover").css({height: height});
    keys_height = parseInt(($(window).height())/2);
    $("section.content:not(:first)").css({height: keys_height*2});
    $("section#home").css({height: keys_height*3+"px"});
    $("section#home").find("div.cell-center").css({marginTop: keys_height+"px"});
    $("section.content:last").css({height: keys_height/10});
    $("div#example-view-container").css({height: keys_height*2-1});
    $("div#example-view").css({height: keys_height*2-1});

    if($("div#example-view").find("pre").length != 0){
        var width = ($(window).width()-900)/2+180;
        $("div#example-view-container").css({width: width});
        $("div#example-view").find("pre").css({height: Math.min(parseInt(content.css("height")), 2*keys_height-60)});
    }
}

$(document).ready(function(){
    
    $("div.example-link").click(function(){
        var content = $("div.example-content#"+$(this).data("target"));
        $("div#example-view").html(content.html());
        $("div.example-link").css({color: "#41B0DB"});
        $(this).css({color: "white"});
        adjustScreen();
    });

    $(document).scroll(function(){
        var position = Math.min($(this).scrollTop(), parseInt(keys_height));
        $("div#cover").css("top", Math.min(-position,0)+"px");
    }); 

    $("div#more").click(function(){ $("body").animate({scrollTop: keys_height}, "slow");     });
    $("div#continue").click(function(){ $("body").animate({scrollTop: 3*keys_height}, "slow");     });

    $("img.sponsor").mouseenter(function(){ $(this).css("-webkit-filter", "grayscale(.0)");  });
    $("img.sponsor").mouseleave(function(){ $(this).css("-webkit-filter", "grayscale(.9)");  });

    adjustScreen(); $(window).resize(adjustScreen);
});

